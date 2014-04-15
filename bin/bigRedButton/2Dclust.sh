#!/bin/bash
#######################
# pairwise structure alignment and clustering
# (C)opyright -- m.smith@garvan.org.au
#
# INPUT: 	output of Bedtools getFasta (in .fasta format, one line per sequence)
#			tested with ">chr1:1000-1100(+)" format in header

#########################################
## CUSTOM PARAMETERS --- MUST BE DEFINED
MAX_CPUS=192
export PROCS=16												# Cpus for clustering and mlocarna
export PATH_TO_SGE_SCRIPTS=${HOME}/apps/bsh/clustering		# scripts used herein
export INPUT_FASTA=$1 										# Input files (see README, if it exists)

### What to run 
RUN_LOCARNA=
RUN_CARNA=
RUN_DOTALIGNER="1"

### load envars
module load gi/ViennaRNA/2.1.3
module load marsmi/dotaligner/27032014
module load marsmi/locarna/1.7.13
# module load marsmi/carna
module load marsmi/newick_utils/1.6
# requires the following R packages: pvclust, snow
module load gi/R/3.0.0

#########################################
# Generate dot plots from fasta file
#>>>>>> Make sure .fasta is one sequence per line
#>>>>>> Ensure name of interest is in fasta file
RNA_FOLD=`which RNAfold`
[[ -z RNA_FOLD ]] && (echo -e "\e[91m[ WARNING ]\e[0m No RNAfold binary in PATH! exiting...\a" ; exit 1 )

## This works on Mac if $1 = `pwd`/file.fa
#FILE_PATH=$INPUT_FASTA 
## This works in linux
FILE_PATH=`readlink -f $INPUT_FASTA`
FILE_NAME=${INPUT_FASTA##*/}
FILE_NAME=${FILE_NAME%.*}
WORK_DIR=${FILE_PATH%/*}

echo -e "\e[93m[ NOTE ]\e[0m Work dir = "$WORK_DIR/$FILE_NAME
echo -e "\e[93m[ NOTE ]\e[0m File name = "$FILE_NAME

## TO DO 
#  check if .ps files already exist (# .ps == # grep ">" fasta )
if [ ! -d  $WORK_DIR/$FILE_NAME ]; 	then	mkdir $WORK_DIR/$FILE_NAME; fi
if [ ! -d  $WORK_DIR/$FILE_NAME/dotplots ]; 
	then
	echo -e "\e[93m[ NOTE ]\e[0m Folding sequences"
	mkdir $WORK_DIR/$FILE_NAME/dotplots
	cd $WORK_DIR/$FILE_NAME/dotplots
	cat $WORK_DIR/$INPUT_FASTA | $RNA_FOLD --noPS --noLP -p
	cd ../
fi
if [[ `ls $WORK_DIR/$FILE_NAME/dotplots/*dp.ps | wc -l` -ne `grep '>' $WORK_DIR/$INPUT_FASTA | wc -l` ]]; 
	then 
	echo -e "\e[91m[ WARNING ]\e[0m inconsistent amount of _dp.ps files --> refolding input"
	echo -e "\e[93m[ NOTE ]\e[0m Folding sequences"

	cd $WORK_DIR/$FILE_NAME/dotplots
	cat $WORK_DIR/$INPUT_FASTA | $RNA_FOLD --noPS --noLP -p 
	cd ../
else
	echo -e "\e[93m[ NOTE ]\e[0m dot plots exist... moving on"
fi

#########################################
# preprocess data for DotAligner input 
echo -e "\e[93m[ NOTE ]\e[0m Parsing dotplots for DotAligner"
for file in $WORK_DIR/$FILE_NAME/dotplots/*dp.ps; do
	if [[ ! -e ${file%ps}pp ]]; then
		fasta_seq=$(grep '\/sequence' -A 1 $file | tail -n 1 | sed 's/\\//g')
		$PATH_TO_SGE_SCRIPTS/getRNAfoldPP.pl $fasta_seq $file > ${file%ps}pp
	fi
done
echo -e "\e[93m[ NOTE ]\e[0m DotAligner input is ready"

#########################################
# outputs all pairwise comparisons from
# a list of file/genes
echo -e "\e[93m[ NOTE ]\e[0m generating list of pairwise comparisons for qsub job arrays"
if [[ ! -e $WORK_DIR/$FILE_NAME/pairwise_comparisons.txt  ]]; then 
	ls $WORK_DIR/$FILE_NAME/dotplots/*dp.ps > $WORK_DIR/$FILE_NAME/file_list.txt
	FILE_LIST=${WORK_DIR}/${FILE_NAME}/file_list.txt
	STRUCTURES=$( wc -l ${FILE_LIST} | awk '{print $1}')
	i=1
	while [ $i -le $STRUCTURES ]; do
 		for (( j = $i ; j <= $STRUCTURES; j++ )); do
	     	echo -n `head -n $i $FILE_LIST | tail -n 1`" "
	     	echo -n `head -n $j $FILE_LIST | tail -n 1`" "
	     	echo $i" "$j
		done
		i=$(($i+1))
	done > $WORK_DIR/$FILE_NAME/pairwise_comparisons.txt 
	# split data files
		LINES=`wc -l ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt | awk '{print $1}'`
		#echo $LINES
		LINES_PER_ARRAY=$(( ${LINES}/${MAX_CPUS}+1  ))
		#echo $LINES_PER_ARRAY
		cd $WORK_DIR/$FILE_NAME
		split -d -a 3 -l $LINES_PER_ARRAY $WORK_DIR/$FILE_NAME/pairwise_comparisons.txt  array_
else
	if [[ ! -e  $WORK_DIR/$FILE_NAME/array_001 ]]; then 
		#########################################
		# split data files
		LINES=`wc -l ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt | awk '{print $1}'`
		#echo $LINES
		LINES_PER_ARRAY=$(( ${LINES}/${MAX_CPUS}+1  ))
		#echo $LINES_PER_ARRAY
		cd $WORK_DIR/$FILE_NAME
		split -d -a 3 -l $LINES_PER_ARRAY $WORK_DIR/$FILE_NAME/pairwise_comparisons.txt  array_
	else 
		echo -e "\e[93m[ NOTE ]\e[0m list of pairwise comparions already exists... moving on "
	fi
fi

#########################################################################
# Launch qsub arrays  
ARRAY_SIZE=$( ls $WORK_DIR/$FILE_NAME/array_* | wc -l )
## DotAligner -- attempts data recovery
if [ ! -z $RUN_DOTALIGNER ]; then 
	if [[ ! -e ${WORK_DIR}/${FILE_NAME}/dotaligner_srtd.out.gz ]]; then 
		echo -e "\e[92m[ QSUB ]\e[0m qsub -cwd -V -N DotAligner -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
				-o ${WORK_DIR}/${FILE_NAME}/DotAligner.o -e ${WORK_DIR}/${FILE_NAME}/DotAligner.e \
				\"${PATH_TO_SGE_SCRIPTS}/dotaligner.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt\""
		DOTALIGNER=$(qsub -cwd -V -N DotAligner -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
				-o ${WORK_DIR}/${FILE_NAME}/DotAligner.o -e ${WORK_DIR}/${FILE_NAME}/DotAligner.e \
				"${PATH_TO_SGE_SCRIPTS}/dotaligner.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt")
	########################################################################
	# POST PROCESSING
		DOTALIGNER_ID=$( echo $DOTALIGNER | cut -d " " -f 3 | cut -d "." -f 1 )
		echo -e "\e[93m[ NOTE ]\e[0m DOTALIGNER job ID: "$DOTALIGNER_ID
		POST_CMD="${PATH_TO_SGE_SCRIPTS}/postprocess.sge ${WORK_DIR}/${FILE_NAME} dotaligner"
		CMD="qsub -hold_jid ${DOTALIGNER_ID} -cwd -V -N DA.postp -pe smp ${PROCS} -b y -j y \
			-o ${WORK_DIR}/${FILE_NAME}/DotAligner.post.o $POST_CMD"
		echo -e "\e[92m[ QSUB ]\e[0m $CMD" && $CMD
	else
		# Delete dotaligner_srtd.out.gz to re-run DotAligner
		echo -e "\e[93m[ NOTE ]\e[0m DotAligner output found! Running post-processing only."
		POST_CMD="${PATH_TO_SGE_SCRIPTS}/postprocess.sge ${WORK_DIR}/${FILE_NAME} dotaligner"
		CMD="qsub -cwd -V -N DA.postp -pe smp ${PROCS} -b y -j y \
			-o ${WORK_DIR}/${FILE_NAME}/DotAligner.post.o $POST_CMD"
		echo -e "\e[92m[ QSUB ]\e[0m $CMD" && $CMD
	fi
fi
## LOCARNA
if [ ! -z $RUN_LOCARNA ]; then 
	echo -e "\e[92m[ QSUB ]\e[0m qsub -cwd -V -N LocaRNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
			-o ${WORK_DIR}/${FILE_NAME}/LocaRNA.o -e ${WORK_DIR}/${FILE_NAME}/LocaRNA.e \ 
			\"${PATH_TO_SGE_SCRIPTS}/locarna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt\""
	LOCARNA=$(qsub -cwd -V -N LocaRNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
			-o ${WORK_DIR}/${FILE_NAME}/LocaRNA.o -e ${WORK_DIR}/${FILE_NAME}/LocaRNA.e \
			"${PATH_TO_SGE_SCRIPTS}/locarna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt")
	# echo $LOCARNA 
	########################################################################
	# POST PROCESSING
	LOCARNA_ID=$( echo $LOCARNA | cut -d " " -f 3 | cut -d "." -f 1 )
	echo -e "\e[93m[ NOTE ]\e[0m LOCARNA job ID: "$LOCARNA_ID
	CMD="qsub -hold_jid ${LOCARNA_ID} -cwd -V -N LCRNA.postp -pe smp 1 -b y -j y \
		-o ${WORK_DIR}/${FILE_NAME}/LocaRNA.post.o \
		\"${PATH_TO_SGE_SCRIPTS}/postprocess.sge ${WORK_DIR}/${FILE_NAME} locarna\""
	echo -e "\e[92m[ QSUB ]\e[0m Post processing: $CMD" && $CMD
fi
## CARNA
if [ ! -z $RUN_CARNA ]; then 
	echo -e "\e[92m[ QSUB ]\e[0m qsub -cwd -V -N CARNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
			-o ${WORK_DIR}/${FILE_NAME}/CARNA.o -e ${WORK_DIR}/${FILE_NAME}/CARNA.e \
			\"${PATH_TO_SGE_SCRIPTS}/carna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt\""
	CARNA=$(qsub -cwd -V -N CARNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
			-o ${WORK_DIR}/${FILE_NAME}/CARNA.o -e ${WORK_DIR}/${FILE_NAME}/CARNA.e \
			"${PATH_TO_SGE_SCRIPTS}/carna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt")
	echo $CARNA
## qsub analyze.sh -W depend=afterokarray:427[]
fi

#########################################
# Post process qsub output 


## qsub -W depend=afterokarray:`echo $LOCARNA | cut -f 2`  postprocess.sge 


## qsub analyze.sh -W depend=afterokarray:427[]


## qsub analyze.sh -W depend=afterokarray:427[]






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
export PATH_TO_SGE_SCRIPTS=${HOME}/apps/bsh/clustering		# scripts used herein
#############################################################
##### Input shoulf be a .fasta file, ideally from "bedtools getFasta" output
##### Please ensure no colons, e.g. ":", are in the fasta header 
##### as this will interfere with .newick tree format 
export INPUT_FASTA=$1 	# Input files (see README, if it exists)
#############################################################
export PROCS=16												# Cpus for clustering and mlocarna
export BOOTSTRAPS=10000										# Number of bootstraps
export ALPHA_STAT=0.99										# Alpha statistic for clustering specificity

### What to run 
RUN_LOCARNA=
RUN_CARNA=
RUN_DOTALIGNER="1"

### load envars
module load gi/ViennaRNA/2.1.3
module load marsmi/dotaligner/27032014
module load marsmi/locarna/1.7.16
# module load marsmi/carna
module load marsmi/newick_utils/1.6
module load stesee/PETfold/prebuilt/2.0
module load gi/R/3.0.0 							# requires the following R packages: pvclust, snow, ape

#########################################
# Generate dot plots from fasta file
#>>>>>> Make sure .fasta is one sequence per line
#>>>>>> Ensure name of interest is in fasta file
RNA_FOLD=`which RNAfold`
[[ -z RNA_FOLD ]] && (echo -e "\e[91m[ ERROR ]\e[0m No RNAfold binary in PATH! exiting...\a" ; exit 1 )

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
	echo -e "\e[93m[ NOTE ]\e[0m Formatting sequence names"
	sed 's/\:/_/' $WORK_DIR/$INPUT_FASTA | sed 's/)//' | sed 's/(/_/' | sed 's/_+/_s/' | sed 's/_-/_r/' | sed 's/-/_/' \
																			 >  $WORK_DIR/$FILE_NAME/formatted_input.fa
	echo -e "\e[93m[ NOTE ]\e[0m Folding sequences"
	mkdir $WORK_DIR/$FILE_NAME/dotplots
	cd $WORK_DIR/$FILE_NAME/dotplots
	cat $WORK_DIR/$FILE_NAME/formatted_input.fa | $RNA_FOLD --noPS --noLP -p
	sleep 2  	# Ensures we don;t get a premature error or warning message from RNAfold not completing on time
	cd ../
fi
if [[ `ls $WORK_DIR/$FILE_NAME/dotplots/*dp.ps | wc -l` -ne `grep '>' $WORK_DIR/$INPUT_FASTA | wc -l` && -e $WORK_DIR/$FILE_NAME/formatted_input.fa ]]; 
	then 
	echo -e "\e[91m[ WARNING ]\e[0m inconsistent amount of _dp.ps files --> refolding input"
	echo -e "\e[93m[ NOTE ]\e[0m Folding sequences"

	cd $WORK_DIR/$FILE_NAME/dotplots
	cat $WORK_DIR/$FILE_NAME/formatted_input.fa | $RNA_FOLD --noPS --noLP -p 
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
# Launch qsub arrays  
ARRAY_SIZE=$( ls ${WORK_DIR}/${FILE_NAME}/array_* | wc -l )
##########################################################################
##################				DOTALIGNER 				##################
if [[ ! -z $RUN_DOTALIGNER ]]; then 
	echo -e "\e[93m[ NOTE ]\e[0m Launching all vs. all pairwise alignments with DotAligner "

## Attempt recovery if alignment was successful
## Delete dotaligner/srtd.out.gz to re-run DotAligner
	if [[ ! -e ${WORK_DIR}/${FILE_NAME}/dotaligner/srtd.out.gz ]]; then 
		# ensure old files are all deleted if bein re-run
		ALN_CMD="${PATH_TO_SGE_SCRIPTS}/dotaligner.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt"
		if [[ -e ${WORK_DIR}/${FILE_NAME}/DotAligner.clust.log ]]; then rm ${WORK_DIR}/${FILE_NAME}/DotAligner.log ; fi
		CMD="qsub -cwd -V -N DotAligner -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y -o ${WORK_DIR}/${FILE_NAME}/DotAligner.log ${ALN_CMD}"
		echo -e "\e[92m[ QSUB ]\e[0m "$CMD && DOTALIGNER_ALN=$( $CMD )
	fi	
##################				CLUSTERING 				##################
## Attempt recovery if clustering was successful
## Delete dotaligner/srtd.newick to re-run clustering
## All clustering will produce a newick tree, but no guarantee there will be clusters 
	if [[ ! -e ${WORK_DIR}/${FILE_NAME}/dotaligner/srtd.newick ]]; then 
		rm -rf ${WORK_DIR}/${FILE_NAME}/dotaligner/* 
		DOTALIGNER_ALN=$( echo ${DOTALIGNER_ALN} | cut -d " " -f 3 | cut -d "." -f 1 )
		echo -e "\e[93m[ NOTE ]\e[0m Clustering output of DOTALIGNER, awaiting completion of job--if required: "$DOTALIGNER_ALN
		echo -e "         This may take a while... you could run this as a background process (ctrl-z; bg)"
		CLUST_CMD="${PATH_TO_SGE_SCRIPTS}/postAlign.sge ${WORK_DIR}/${FILE_NAME} dotaligner"
		if [[ -e ${WORK_DIR}/${FILE_NAME}/DotAligner.clust.log ]]; then rm ${WORK_DIR}/${FILE_NAME}/DotAligner.post.log ; fi
		CMD="qsub -sync y -hold_jid ${DOTALIGNER_ALN} -cwd -V -N DA.clust -pe smp ${PROCS} -b y -j y \
			-o ${WORK_DIR}/${FILE_NAME}/DotAligner.clust.log $CLUST_CMD"
		echo -e "\e[92m[ QSUB ]\e[0m $CMD" && DOTALIGNER_CLUST=$( $CMD )
	fi

##################				POSTPROCESSING 				##################
## Post processing of any identified clusters. 
#### ${1}/${2}/srtd_a${ALPHA_STAT}_clusters/cluster.X/cluster.X.fa 
#	if [[ ! -e ${WORK_DIR}/${FILE_NAME}/clusters/cluster.1/cluster.1.fa ]]; then
	DOTALIGNER_CLUST=$( echo ${DOTALIGNER_CLUST} | cut -d " " -f 3 | cut -d "." -f 1 )		# Get job ID of existing process
	ARRAY_SIZE=$( ls ${WORK_DIR}/${FILE_NAME}/dotaligner/clusters/cluster.*.fa | wc -l ) 						# Get amount of clusters to process
	if [[ ARRAY_SIZE -eq 0 ]]; then 
		echo -e "\e[91m[ WARNING ]\e[0m No clusters found in ${WORK_DIR}/${FILE_NAME}/dotaligner/clusters/"
	fi
	POST_CMD="${PATH_TO_SGE_SCRIPTS}/postClust.sge ${WORK_DIR}/${FILE_NAME} dotaligner"		# Setup the command and args
	if [[ -e ${WORK_DIR}/${FILE_NAME}/DotAligner.post.log ]]; then rm ${WORK_DIR}/${FILE_NAME}/DotAligner.post.log ; fi
	CMD="qsub -hold_jid ${DOTALIGNER_CLUST} -cwd -V -N DA.postp -pe smp ${PROCS} -t 1-${ARRAY_SIZE} -b y -j y -o ${WORK_DIR}/${FILE_NAME}/DotAligner.post.log ${POST_CMD}"						# Setup SGE command
	
	echo -e "\e[92m[ QSUB ]\e[0m $CMD" && DOTALIGNER_POST=$( $CMD )							# Print command and execute
fi


### CLEAN UP NOTES
# rm ${WORK_DIR}/${FILE_NAME}/array_*.locarna
# rm ${WORK_DIR}/${FILE_NAME}/array_*

# ## LOCARNA
# if [ ! -z $RUN_LOCARNA ]; then 
# 	echo -e "\e[92m[ QSUB ]\e[0m qsub -cwd -V -N LocaRNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
# 			-o ${WORK_DIR}/${FILE_NAME}/LocaRNA.o -e ${WORK_DIR}/${FILE_NAME}/LocaRNA.e \ 
# 			\"${PATH_TO_SGE_SCRIPTS}/locarna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt\""
# 	LOCARNA=$(qsub -cwd -V -N LocaRNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
# 			-o ${WORK_DIR}/${FILE_NAME}/LocaRNA.o -e ${WORK_DIR}/${FILE_NAME}/LocaRNA.e \
# 			"${PATH_TO_SGE_SCRIPTS}/locarna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt")
# 	# echo $LOCARNA 
# 	########################################################################
# 	# POST PROCESSING
# 	LOCARNA_ID=$( echo $LOCARNA | cut -d " " -f 3 | cut -d "." -f 1 )
# 	echo -e "\e[93m[ NOTE ]\e[0m LOCARNA job ID: "$LOCARNA_ID
# 	CMD="qsub -hold_jid ${LOCARNA_ID} -cwd -V -N LCRNA.postp -pe smp 1 -b y -j y \
# 		-o ${WORK_DIR}/${FILE_NAME}/LocaRNA.post.o \
# 		\"${PATH_TO_SGE_SCRIPTS}/postAlign.sge ${WORK_DIR}/${FILE_NAME} locarna\""
# 	echo -e "\e[92m[ QSUB ]\e[0m Post processing: $CMD" && $CMD
# fi
# ## CARNA
# if [ ! -z $RUN_CARNA ]; then 
# 	echo -e "\e[92m[ QSUB ]\e[0m qsub -cwd -V -N CARNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
# 			-o ${WORK_DIR}/${FILE_NAME}/CARNA.o -e ${WORK_DIR}/${FILE_NAME}/CARNA.e \
# 			\"${PATH_TO_SGE_SCRIPTS}/carna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt\""
# 	CARNA=$(qsub -cwd -V -N CARNA -pe smp 1 -t 1-${ARRAY_SIZE} -b y -j y \
# 			-o ${WORK_DIR}/${FILE_NAME}/CARNA.o -e ${WORK_DIR}/${FILE_NAME}/CARNA.e \
# 			"${PATH_TO_SGE_SCRIPTS}/carna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt")
# 	echo $CARNA
# ##
#  qsub analyze.sh -W depend=afterokarray:427[]
# fi

# #########################################
# # Post process qsub output 
# ## qsub -W depend=afterokarray:`echo $LOCARNA | cut -f 2`  postAlign.sge 
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
export PATH_TO_SGE_SCRIPTS=${HOME}/DotAligner/bin/bigRedButton		# scripts used herein
export GENOME=/share/ClusterShare/biodata/contrib/genomeIndices_garvan/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa      # fasta file of reference genome (here hg19)
export PROCS=8												# Cpus for clustering and mlocarna
export BOOTSTRAPS=10000										# Number of bootstraps
export ALPHA_STAT=0.99										# Alpha statistic for clustering specificity
export BETA_STAT=0.0
export SPAN=150

### What to run 
RUN_LOCARNA=
RUN_CARNA=
RUN_DOTALIGNER="1"

#########################################
# DotAligner parameters
KAPPA=0.3
ALPHA=-0.05 #-0.1
BETA=-0.05 #-0.001
RADIUS=1 #5
THETA=0.5
DELTANULL=0 #0.5
SEEDLEN=0
MAXSHIFT=10 #20
SEQALN="" #"--seqaln"                        #set to "--seqaln" or ""
PRECISION=4 #2
PNULL=0.0005 #0.01

#########################################
# usage
usage() {
	echo Program: "2Dclust.sh (pairwise structure alignment and clustering)"
	echo Author: Stefan Seemann, Martin Smith
	echo Version: 1.0
	echo Contact: m.smith@garvan.org.au
	echo "Usage: 2Dclust.sh -b <file> [OPTIONS]"
	echo "         or"
	echo "       2Dclust.sh -f <file> [OPTIONS]"
	echo "Options:"
	echo " -b <file>    .. input file of signals (BED format)"	
	echo "                 -> fasta file is made for you of sequences representing local windows"
	echo "                    of highest structure mass that overlap the signals"
	echo " -f <file>    .. input file of sequences that are output of Bedtools getFasta (FASTA format, one line per sequence)"
	echo "                 -> these sequences are direct input for clustering"
	echo " -w <length>  .. length of local windows considered for RNA structure clustering"
	echo " -h           .. help"
	echo
	exit 0
}

#########################################
# parse options
while getopts b:f:w:h ARG; do
	case "$ARG" in
		#############################################################
		##### Input shoulf be a .fasta file, ideally from "bedtools getFasta" output
		##### Please ensure no colons, e.g. ":", are in the fasta header 
		##### as this will interfere with .newick tree format 
		b) export INPUT_BED=$OPTARG;;
		f) export INPUT_FASTA=$OPTARG;; 
		w) SPAN=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
[ "$INPUT_BED" -a ! -f "$INPUT_BED" ] && usage;
[ "$INPUT_FASTA" -a ! -f "$INPUT_FASTA" ] && usage;
[ "$HELP" ] && usage;
[ "$INPUT_BED" -a "$INPUT_FASTA" ] && usage;
[ ! "$INPUT_BED" -a ! "$INPUT_FASTA" ] && usage;


#########################################
# load envars 
# N.B. Make sure all required binaires are in yout $PATH if you don't use rocks-modules
module load gi/ViennaRNA/2.1.3
#module load marsmi/dotaligner/27032014
module load marsmi/dotaligner/0.2
module load marsmi/locarna/1.7.16
module load marsmi/newick_utils/1.6
module load stesee/PETfold/prebuilt/2.0
module load stesee/RNAbound/prebuilt/1.1
module load gi/R/3.0.0 							# requires the following R packages: pvclust, snow, ape
module load gi/bedtools/2.19.1
# module load marsmi/carna



#########################################
# read command line arguments
## This works on Mac if $1 = `pwd`/file.fa
#FILE_PATH=$INPUT_FASTA 
## This works in linux
#FILE_PATH=`readlink -f $INPUT_FASTA`
if [ -f "$INPUT_FASTA" ];
then
	INPUT=$INPUT_FASTA
	if [ `echo $INPUT_FASTA | awk '$1!~/^\//'` ]; then INPUT_FASTA=$PWD/$INPUT_FASTA; fi
else
	INPUT=$INPUT_BED
	if [ `echo $INPUT_BED | awk '$1!~/^\//'` ]; then INPUT_BED=$PWD/$INPUT_BED; fi
fi
TEMP_NAME=${INPUT##*/}
TEMP_NAME=${TEMP_NAME%.*}
WORK_DIR=$PWD
if [ -f "$INPUT_FASTA" ]; then FILE_NAME=$TEMP_NAME; else FILE_NAME=${TEMP_NAME}_bound_signal_span${SPAN}; fi

echo -e "\e[93m[ NOTE ]\e[0m Work dir = "$WORK_DIR/$FILE_NAME
echo -e "\e[93m[ NOTE ]\e[0m File name = "$FILE_NAME
if [[ ! -d  $WORK_DIR/$FILE_NAME ]]; 	then	mkdir $WORK_DIR/$FILE_NAME; fi
if [[ -f "$INPUT_FASTA" && ! -e $WORK_DIR/$FILE_NAME/${FILE_NAME}.fasta ]]
	then ln -s $INPUT_FASTA $WORK_DIR/$FILE_NAME/${FILE_NAME}.fasta
fi


#########################################
# find local window of highest structure mass around input signals
if [ -f "$INPUT_BED" ];
then 
	if [ ! -d  ${WORK_DIR}/${FILE_NAME}/preprocesseddata ]; 
	then
		mkdir ${WORK_DIR}/${FILE_NAME}/preprocesseddata
		cd ${WORK_DIR}/${FILE_NAME}/preprocesseddata
	
		#merge and extend
		echo -e "\e[93m{ NOTE ]\e[0m Merge input signals of less than 50nt distance"
		bedtools merge -scores max -s -d 50 -nms -i ${INPUT_BED} > ${TEMP_NAME}_sorted_merged.bed  
		WINDOW=$((SPAN+SPAN))
       		echo -e "\e[93m{ NOTE ]\e[0m Create windows of $WINDOW nucleotides centered by merged input signals"
 		awk -v span=$SPAN 'BEGIN{OFS="\t"}{m=$2+int(($3-$2)/2); print $1,m-span,m+span,$4,$5,$6}' ${TEMP_NAME}_sorted_merged.bed > ${TEMP_NAME}_sorted_extended.bed

		#get sequences
       		echo -e "\e[93m{ NOTE ]\e[0m Get fasta file"
		bedtools getfasta -fi $GENOME -bed ${TEMP_NAME}_sorted_extended.bed -s -fo ${TEMP_NAME}_sorted_extended.fasta
	fi

	if [ ! -d  ${WORK_DIR}/${FILE_NAME}/rnabound ]; 
	then
		RNA_BOUND=`which RNAbound`
		[[ -z RNA_BOUND ]] && (echo -e "\e[91m[ ERROR ]\e[0m No RNAbound binary in PATH! exiting...\a" ; exit 1 )

		WIN=$((SPAN+50))
		echo -e "\e[93m{ NOTE ]\e[0m Find local windows of length $SPAN that have highest structure mass"
		mkdir ${WORK_DIR}/${FILE_NAME}/rnabound
		cd ${WORK_DIR}/${FILE_NAME}/rnabound
	
	  	#create dotplot using RNAplfold
	   	RNAplfold -W $WIN -L $SPAN -c 0.0005 < ${WORK_DIR}/${FILE_NAME}/preprocesseddata/${TEMP_NAME}_sorted_extended.fasta
	
		for DOTPLOT in *_dp.ps;
		do
			PNAME=${DOTPLOT%_*}
	
		      	#convert the dotplot into matrix
		 	dot2matrix.pl $DOTPLOT > ${PNAME}.bpp
	
			#execute RNAbound using the matrix file
			#RNAbound -M ${PNAME}.bpp -g1 -w $SPAN | subsetrnabound.pl - 0.5 > ${PNAME}.bound
			RNAbound -M ${PNAME}.bpp -g1 -w $SPAN | head -1 > ${PNAME}.bound
			awk -v name=$PNAME 'BEGIN{OFS="\t"}{split($1,a,"-"); split(name,b,":"); split(b[2],c,"("); split(c[1],d,"-"); sub(")","",c[2]); print b[1],a[1]-1+d[1],a[2]+d[1],name,"0",c[2]}' ${PNAME}.bound >> ${TEMP_NAME}.bound_span${SPAN}.bed
		done
	
		#check that window overlaps input signal (CHECK!!!)
		bedtools intersect -wa -a ${TEMP_NAME}.bound_span${SPAN}.bed -b ${INPUT_BED} | sort -u > ${FILE_NAME}.bed
	
		#get sequence of local structure
		bedtools getfasta -fi $GENOME -bed ${FILE_NAME}.bed -name -fo ${FILE_NAME}.fasta
		ln -s ${WORK_DIR}/${FILE_NAME}/rnabound/${FILE_NAME}.fasta $WORK_DIR/$FILE_NAME/${FILE_NAME}.fasta
	
		rm -f *.bpp *_dp.ps *.bound.bed *.bound
	fi
fi

#########################################
# Generate dot plots from fasta file
#>>>>>> Make sure .fasta is one sequence per line
#>>>>>> Ensure name of interest is in fasta file
## TO DO 
#  check if .ps files already exist (# .ps == # grep ">" fasta )
RNA_FOLD=`which RNAfold`
[[ -z RNA_FOLD ]] && (echo -e "\e[91m[ ERROR ]\e[0m No RNAfold binary in PATH! exiting...\a" ; exit 1 )

if [ ! -d  $WORK_DIR/$FILE_NAME/dotplots ]; 
then
	echo -e "\e[93m[ NOTE ]\e[0m Formatting sequence names"
	sed 's/\:/_/' ${WORK_DIR}/${FILE_NAME}/${FILE_NAME}.fasta | sed 's/)//' | sed 's/(/_/' | sed 's/_+/_s/' | sed 's/_-/_r/' | sed 's/-/_/' \
																			 >  $WORK_DIR/$FILE_NAME/formatted_input.fa
	echo -e "\e[93m[ NOTE ]\e[0m Folding sequences"
	mkdir $WORK_DIR/$FILE_NAME/dotplots
	cd $WORK_DIR/$FILE_NAME/dotplots
	cat $WORK_DIR/$FILE_NAME/formatted_input.fa | $RNA_FOLD --noPS --noLP -p
	sleep 2  	# Ensures we don;t get a premature error or warning message from RNAfold not completing on time
	cd ../
fi
if [[ `ls $WORK_DIR/$FILE_NAME/dotplots/*dp.ps | wc -l` -ne `grep '>' $WORK_DIR/$FILE_NAME/${FILE_NAME}.fasta | wc -l` && -e $WORK_DIR/$FILE_NAME/formatted_input.fa ]]; 
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
			
			## This messes up downstream processing in postAlign.sge
			## It produces :
			## 	Filename_1 Filename_2 Filename_1 Filename_2
			## instead of :
			##	Filename_1 Filename_2 1 2
			## vvvvvvvvvvvvvvvvvvvvvvvvvvvvv

			#FILE1=`head -n $i $FILE_LIST | tail -n 1`
			#FILE2=`head -n $j $FILE_LIST | tail -n 1`
			#TEMP=${FILE1%___*}
			#IDX1=${TEMP##*___}
			#TEMP=${FILE2%___*}
			#IDX2=${TEMP##*___}
			#echo $FILE1" "$FILE2" "$IDX1" "$IDX2  
			
			## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
		LINES_PER_ARRAY=$(( ${LINES}/${MAX_CPUS}+1  )) # ensures less array jobs are run than max_slots; balanced IO/CPU performance 
		cd $WORK_DIR/$FILE_NAME
		split -d -a 3 -l $LINES_PER_ARRAY $WORK_DIR/$FILE_NAME/pairwise_comparisons.txt  array_
	else 
		echo -e "\e[93m[ NOTE ]\e[0m list of pairwise comparions already exists... moving on "
	fi
fi
# Launch qsub arrays  
ARRAY_SIZE=$( ls ${WORK_DIR}/${FILE_NAME}/array_* | wc -l )
echo -e "\e[93m[ NOTE ]\e[0m Launching "$ARRAY_SIZE" arrays of "$LINES_PER_ARRAY" jobs"

##########################################################################
##################				DOTALIGNER 				##################
if [[ ! -z $RUN_DOTALIGNER ]]; then 
## Attempt recovery if alignment was successful
## Delete dotaligned.checkpoint to re-run DotAligner
	if [[ ! -e ${WORK_DIR}/${FILE_NAME}/dotaligned.checkpoint ]]; then 
		echo -e "\e[93m[ NOTE ]\e[0m Launching all vs. all pairwise alignments with DotAligner "
		# ensure old files are all deleted if bein re-run
		#ALN_CMD="${PATH_TO_SGE_SCRIPTS}/dotaligner.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt"
		ALN_CMD="${PATH_TO_SGE_SCRIPTS}/dotaligner.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt $KAPPA $ALPHA $BETA $RADIUS $THETA $DELTANULL $SEEDLEN $MAXSHIFT $PRECISION $PNULL $SEQALN"
		if [[ -e ${WORK_DIR}/${FILE_NAME}/DotAligner.clust.log ]]; then rm ${WORK_DIR}/${FILE_NAME}/DotAligner.log ; fi
		echo -e "\e[93m[ NOTE ]\e[0m DotAligner parameters:"
		echo -e "     KAPPA = "$KAPPA"\nALPHA = "$ALPHA"\nBETA = "$BETA"\nRADIUS = "$RADIUS"\nTHETA = "$THETA"\nDELTANULL = "$DELTANULL"\nSEEDLEN = "$SEEDLEN"\nMAXSHIFT = "$MAXSHIFT"\nPRECISION = "$PRECISION"\nPNULL = "$PNULL"\nSEQALN = "$SEQALN > ${WORK_DIR}/${FILE_NAME}/DotAligner.log
		CMD="qsub -cwd -V -N DotAligner -pe smp 1 -l h_vmem=1G -t 1-${ARRAY_SIZE} -b y -j y -o ${WORK_DIR}/${FILE_NAME}/DotAligner.log ${ALN_CMD} && touch ${WORK_DIR}/${FILE_NAME}/dotaligned.checkpoint"
		echo -e "\e[92m[ QSUB ]\e[0m "$CMD && DOTALIGNER_ALN=$( $CMD )
	else
		echo -e "\e[93m[ NOTE ]\e[0m Parwise alignments already exist! Moving on... "
		echo -e "         If you want to re-run this step, please delete checkpoint file: "${WORK_DIR}/${FILE_NAME}/dotaligned.checkpoint 
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
	if [[ ! -e ${WORK_DIR}/${FILE_NAME}/clusters/mlocarna.cluster.1/results/result.aln ]]; then
		DOTALIGNER_CLUST=$( echo ${DOTALIGNER_CLUST} | cut -d " " -f 3 | cut -d "." -f 1 )		# Get job ID of existing process
		ARRAY_SIZE=$( ls ${WORK_DIR}/${FILE_NAME}/dotaligner/clusters/cluster.*.fa | wc -l ) 						# Get amount of clusters to process
		if [[ ARRAY_SIZE -eq 0 ]]; then 
			echo -e "\e[91m[ WARNING ]\e[0m No clusters found in ${WORK_DIR}/${FILE_NAME}/dotaligner/clusters/"
		else
			echo -e "\e[93m[ NOTE ]\e[0m ${ARRAY_SIZE} clusters have been found!"
		fi
		POST_CMD="${PATH_TO_SGE_SCRIPTS}/postClust.sge ${WORK_DIR}/${FILE_NAME} dotaligner"		# Setup the command and args
		if [[ -e ${WORK_DIR}/${FILE_NAME}/DotAligner.post.log ]]; then rm ${WORK_DIR}/${FILE_NAME}/DotAligner.post.log ; fi
		CMD="qsub -hold_jid ${DOTALIGNER_CLUST} -cwd -V -N DA.postp -pe smp ${PROCS} -t 1-${ARRAY_SIZE} -b y -j y -o ${WORK_DIR}/${FILE_NAME}/DotAligner.post.log ${POST_CMD}"		# Setup SGE command
		echo -e "\e[92m[ QSUB ]\e[0m $CMD" && DOTALIGNER_POST=$( $CMD )							# Print command and execute
	fi
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

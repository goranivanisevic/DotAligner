#!/bin/bash
#######################
# pairwise structure alignment and clustering
# (C)opyright -- m.smith@garvan.org.au
#
MAX_CPUS=192
PATH_TO_SGE_SCRIPTS=${HOME}/apps/bsh/clustering
INPUT_FASTA=$1
#########################################
# Generate dot plots from fasta file
#>>>>>> Make sure .fasta is one sequence per line
#>>>>>> Ensure name of interest is in fasta file

module load gi/ViennaRNA/2.1.3
RNA_FOLD=`which RNAfold`

## This works on Mac if $1 = `pwd`/file.fa
#FILE_PATH=$INPUT_FASTA 
## This works in linux
FILE_PATH=`readlink -f $INPUT_FASTA`
FILE_NAME=${INPUT_FASTA##*/}
FILE_NAME=${INPUT_FASTA%.*}
WORK_DIR=${FILE_PATH%/*}

echo "Work dir = "$WORK_DIR/$FILE_NAME
echo "File name = "$FILE_NAME

## TO DO 
#  check if .ps files already exist (# .ps == # grep ">" fasta )
if [ ! -d  $WORK_DIR/$FILE_NAME ]; 	then	mkdir $WORK_DIR/$FILE_NAME; fi
if [ ! -d  $WORK_DIR/$FILE_NAME/dotplots ]; 
	then
	echo "[ NOTE ] Folding sequences"
	mkdir $WORK_DIR/$FILE_NAME/dotplots
	cd $WORK_DIR/$FILE_NAME/dotplots
	cat $WORK_DIR/$INPUT_FASTA | $RNA_FOLD --noPS --noLP -p
	cd ../
fi
if [[ `ls $WORK_DIR/$FILE_NAME/dotplots/*dp.ps | wc -l` -ne `grep '>' $WORK_DIR/$INPUT_FASTA | wc -l` ]]; 
	then 
	echo "[ WARNING ] inconsistent amount of _dp.ps files --> refolding input"
	echo "[ NOTE ] Folding sequences"

	cd $WORK_DIR/$FILE_NAME/dotplots
	cat $WORK_DIR/$INPUT_FASTA | $RNA_FOLD --noPS --noLP -p
	cd ../
fi

#########################################
# outputs all pairwise comparisons from
# a list of file/genes
ls $WORK_DIR/$FILE_NAME/dotplots/*dp.ps > $WORK_DIR/$FILE_NAME/file_list.txt
FILE_LIST=${WORK_DIR}/${FILE_NAME}/file_list.txt
STRUCTURES=$( wc -l ${FILE_LIST} | awk '{print $1}')
i=1
while [ $i -lt $STRUCTURES ]; do
 	for (( j = $i+1 ; j <= $STRUCTURES; j++ )); do
     	echo -n `head -n $i $FILE_LIST | tail -n 1`" "
     	echo -n `head -n $j $FILE_LIST | tail -n 1`" "
     	echo $i" "$j
	done
	i=$(($i+1))
done > $WORK_DIR/$FILE_NAME/pairwise_comparisons.txt 

#########################################
# split data files
# 
LINES=`wc -l ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt | awk '{print $1}'`
#echo $LINES
LINES_PER_ARRAY=$(( ${LINES}/${MAX_CPUS}+1  ))
#echo $LINES_PER_ARRAY
cd $WORK_DIR/$FILE_NAME
split -d -a 3 -l $LINES_PER_ARRAY $WORK_DIR/$FILE_NAME/pairwise_comparisons.txt  array_


#########################################
# Launch qsub array  
# 
## LOCARNA
echo "qsub -cwd -V -N LocaRNA -pe smp 1 -t 1-192 -b y \"${PATH_TO_SGE_SCRIPTS}/locarna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt\""
## CARNA
echo "qsub -cwd -V -N CARNA -pe smp 1 -t 1-192 -b y \"${PATH_TO_SGE_SCRIPTS}/carna.sge ${WORK_DIR}/${FILE_NAME}/pairwise_comparisons.txt\""


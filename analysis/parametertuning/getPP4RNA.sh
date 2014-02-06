#!/bin/bash
#
# generate seperate fasta file and pp file (basepair probabilities calculated by 'RNAfold -p -d2 --noLP')
# from an input fasta file listing many sequences (sequences span only one line per identifier)
#

PDOTALIGNER='/media/data/phd/projects/DotAligner'
#RNAINPUT=${PDOTALIGNER}/analysis/rfam/selected_PIDs/seed_56_95.fasta
#RNAINPUT=${PDOTALIGNER}/analysis/benchmarking_seqs/sample_pairwise_aligns.fasta
RNAINPUT=$1

for ID in `cat ${RNAINPUT} | awk '/^>/'`;
do
  grep -A 1 ${ID} ${RNAINPUT} | head -2 | awk '{if(/^>/){split($1,a,"/"); gsub(/\./,"_",a[1]); print a[1]}else{print $0}}' > tmp.fasta;
  SEQNAME=`head -1 tmp.fasta | sed 's/^>//'`;
  mv tmp.fasta ${SEQNAME}.fasta;
  ${PDOTALIGNER}/bin/getRNAfoldPP.pl ${SEQNAME}.fasta > ${SEQNAME}.pp;
done


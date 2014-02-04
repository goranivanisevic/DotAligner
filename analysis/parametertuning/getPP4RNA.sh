PDOTALIGNER='/media/data/phd/projects/DotAligner'
RNAINPUT=${PDOTALIGNER}/analysis/rfam/selected_PIDs/seed_56_95.fasta

for ID in `cat ${RNAINPUT} | awk '/^>/'`;
do
  grep -A 1 ${ID} ${RNAINPUT} | head -2 | awk '{if(/^>/){split($1,a,"/"); gsub(/\./,"_",a[1]); print a[1]}else{print $0}}' > tmp.fasta;
  SEQNAME=`head -1 tmp.fasta | sed 's/^>//'`;
  mv tmp.fasta ${SEQNAME}.fasta;
  ${PDOTALIGNER}/bin/getRNAfoldPP.pl ${SEQNAME}.fasta > ${SEQNAME}.pp;
done


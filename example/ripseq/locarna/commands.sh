HDOTALIGNER=/home/stesee/DotAligner
export PROG=locarna
export DATA=DNMT1

#merge results
gunzip -c locarna.out.gz | awk 'BEGIN{OFS=","}{if($1=="ALIGNING"){split($2,a,"/");split($3,b,"/");sub(/_dp.ps/,"",a[length(a)]);sub(/_dp.ps/,"",b[length(a)]);s1=a[length(a)]; s2=b[length(b)]};if($1=="Score:"){sc=$2; print s1,s2,sc}}' | gzip > ${PROG}.${DATA}.csv.gz
 
#pvclust
gunzip -c ${PROG}.${DATA}.csv.gz | awk '{gsub(/,/,"\t");print $0}' | awk '{print $1"\n"$2}' | sort -u > ${DATA}.identifier
LEN=`cat ${DATA}.identifier | wc -l`
gunzip -c ${PROG}.${DATA}.csv.gz | awk '{gsub(/,/,"\t");print $0}' | awk '{print $1,$2,$3"\n"$2,$1,$3}' > tmp.${DATA}.${PROG}.tsv
gunzip -c ${PROG}.${DATA}.csv.gz | awk '{gsub(/,/,"\t");print $0}' | awk '{print $1"\n"$2}' | sort -u | awk '{print $1,$1,20000}' >> tmp.${DATA}.${PROG}.tsv
cat tmp.${DATA}.${PROG}.tsv | sort -k1,2 | awk -v len=$LEN 'BEGIN{ORS=""}{print $3," "; if(NR%len == 0){print "\n"}}' > ${DATA}.${PROG}.2x2
#submit job
qsub -q all.q -N pc4${PROG} job_launch_pvclust.sge

#get members of clusters
cat ${DATA}.${PROG}.clusters | awk '{if(/^\[\[/){gsub(/\[/,"",$1);gsub(/\]/,"",$1);c=$1}else{if(/\[\s*/){for(i=2;i<=NF;i++){gsub(/"/,"",$i);print c,$i}}}}' > ${DATA}.${PROG}.clusters.list
#get non-clustered families as FN
merge.pl <(cut -d" " -f2 ${DATA}.${PROG}.clusters.list | sort | uniq -c) <(ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/${DATA}/*.fa | awk '{split($1,a,/\//);split(a[length(a)],b,/_/);print b[1]}' | sort | uniq -c ) 2,2 | awk '{print ($1-$3),$2}' > ${DATA}.${PROG}.noclusters.list
${HDOTALIGNER}/bin/clusterstatistics.pl ${DATA}.${PROG}.clusters.list ${DATA}.${PROG}.noclusters.list

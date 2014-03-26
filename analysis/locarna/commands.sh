HDOTALIGNER=/home/stesee/DotAligner
export PROG=locarna
export DATA=seed_10_55.new
#export DATA=seed_56_95.new

egrep -A 2 -e '^(\-{10})$' seed_56_95_new_merged.out |  while read line 
do if [[ `echo $line | cut -c 1-3` == "ALI" ]]
	then 
		echo $line | sed 's/\/home\/marsmi\/data\/clustering\/seed_56_95_new\/dotplots\///g' | \
					 awk '{printf $2"\t"; printf $3"\t"; printf $4"\t"; printf $5}' 
elif [[ `echo $line | cut -c 1-3` == "Sco" ]]
	then 
		echo $line | awk '{print "\t"$2}'  
fi 
done > seed_56_95_new_merged.tsv
#cat seed_56_95_new_merged.out | grep -C 1 ALIG | while read line ; do if [[ `echo $line | cut -c 1-3` == "ALI" ]]; then echo $line | sed 's/\/home\/marsmi\/data\/clustering\ntf $2"\t"; printf $3"\t"; printf $4"\t"; printf $5}' ; elif [[ `echo $line | cut -c 1-3` == "Sco" ]]; then echo $line | awk '{print "\t"$2}'  ; fi ; done > seed_56_95_new_merged.tsv

#pvclust
cat ${DATA}_merged.tsv | awk '{print $1"\n"$2}' | sort -u | awk '{split($1,a,"_");print a[1]}' > ${DATA}.identifier
LEN=`cat ${DATA}.identifier | wc -l`
cat ${DATA}_merged.tsv | awk '{print $1,$2,$5"\n"$2,$1,$5}' > tmp.${DATA}.merged.tsv
cat ${DATA}_merged.tsv | awk '{print $1"\n"$2}' | sort -u | awk '{print $1,$1,20000}' >> tmp.${DATA}.merged.tsv
cat tmp.${DATA}.merged.tsv | sort -k1,2 | awk -v len=$LEN 'BEGIN{ORS=""}{print $3," "; if(NR%len == 0){print "\n"}}' > ${DATA}.${PROG}.2x2
#submit job
qsub -q all.q -N pc4${PROG} job_launch_pvclust.sge

#get members of clusters
cat ${DATA}.${PROG}.clusters | awk '{if(/^\[\[/){gsub(/\[/,"",$1);gsub(/\]/,"",$1);c=$1}else{if(/\[\s*/){for(i=2;i<=NF;i++){gsub(/"/,"",$i);print c,$i}}}}' > ${DATA}.${PROG}.clusters.list
#get non-clustered families as FN
merge.pl <(cut -d" " -f2 ${DATA}.${PROG}.clusters.list | sort | uniq -c) <(ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/${DATA}/*.fa | awk '{split($1,a,/\//);split(a[length(a)],b,/_/);print b[1]}' | sort | uniq -c ) 2,2 | awk '{print ($1-$3),$2}' > ${DATA}.${PROG}.noclusters.list
${HDOTALIGNER}/bin/clusterstatistics.pl ${DATA}.${PROG}.clusters.list ${DATA}.${PROG}.noclusters.list

DATA=seed_10_55.new
RF01699:	TP = 4	FP = 0	TN = 134	FN = 13	SP = 1.0000	SN = 0.2353
RF02001:	TP = 1	FP = 2	TN = 129	FN = 19	SP = 0.9847	SN = 0.0500
RF00015:	TP = 11	FP = 15	TN = 122	FN = 3	SP = 0.8905	SN = 0.7857
RF02003:	TP = 9	FP = 0	TN = 142	FN = 0	SP = 1.0000	SN = 1.0000
RF00059:	TP = 13	FP = 4	TN = 127	FN = 7	SP = 0.9695	SN = 0.6500
RF00379:	TP = 6	FP = 0	TN = 139	FN = 6	SP = 1.0000	SN = 0.5000
RF01055:	TP = 14	FP = 12	TN = 120	FN = 5	SP = 0.9091	SN = 0.7368
RF01794:	TP = 3	FP = 0	TN = 141	FN = 7	SP = 1.0000	SN = 0.3000
RF00557:	TP = 2	FP = 0	TN = 131	FN = 18	SP = 1.0000	SN = 0.1000
RF00162:	TP = 3	FP = 14	TN = 133	FN = 1	SP = 0.9048	SN = 0.7500
RF00380:	TP = 5	FP = 0	TN = 145	FN = 1	SP = 1.0000	SN = 0.8333
SP = 0.9690
SN = 0.5401

DATA=seed_56_95.new
RF01699:	TP = 7	FP = 59	TN = 88	FN = 0	SP = 0.5986	SN = 1.0000
RF02001:	TP = 9	FP = 57	TN = 87	FN = 1	SP = 0.6042	SN = 0.9000
RF00015:	TP = 7	FP = 0	TN = 134	FN = 13	SP = 1.0000	SN = 0.3500
RF01685:	TP = 19	FP = 0	TN = 134	FN = 1	SP = 1.0000	SN = 0.9500
RF00379:	TP = 7	FP = 0	TN = 134	FN = 13	SP = 1.0000	SN = 0.3500
RF01055:	TP = 14	FP = 0	TN = 140	FN = 0	SP = 1.0000	SN = 1.0000
RF00557:	TP = 2	FP = 0	TN = 151	FN = 1	SP = 1.0000	SN = 0.6667
RF00162:	TP = 3	FP = 63	TN = 71	FN = 17	SP = 0.5299	SN = 0.1500
RF00380:	TP = 20	FP = 46	TN = 88	FN = 0	SP = 0.6567	SN = 1.0000
RF00013:	TP = 20	FP = 46	TN = 88	FN = 0	SP = 0.6567	SN = 1.0000
SP = 0.8046
SN = 0.7367


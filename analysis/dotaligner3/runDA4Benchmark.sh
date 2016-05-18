BIN=/home/stesee/bin
HDOTALIGNER=/home/stesee/DotAligner
PROG=dotaligner3
EXE=${HDOTALIGNER}/bin/DotAligner

DATA=$1
#DATA=seed_10_55
#DATA=seed_56_95

rm -f ${DATA}.${PROG};
for FAMSEQ1 in `ls ${HDOTALIGNER}/analysis/${PROG}/data/${DATA}/*.pp`;
do
  for FAMSEQ2 in `ls ${HDOTALIGNER}/analysis/${PROG}/data/${DATA}/*.pp`;
  do
    echo $FAMSEQ1 > tmp.famseq
    echo $FAMSEQ2 >> tmp.famseq
    if [[ `sort tmp.famseq | head -1` == $FAMSEQ1 ]];
    then

      echo "echo $FAMSEQ1 $FAMSEQ2 >> ${DATA}.${PROG};"
      echo "{ time ${EXE} -d $FAMSEQ1 -d $FAMSEQ2 -k 0.3 -a -0.05 -b -0.05 -r 5 -t 0.5 -l 0.5 -s 0 -m 15; } &>> ${DATA}.${PROG};"

    fi;
  done;
done;


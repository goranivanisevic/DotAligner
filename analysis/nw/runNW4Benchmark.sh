BIN=/home/stesee/bin
HDOTALIGNER=/home/stesee/DotAligner
PROG=nw
EXE=${HDOTALIGNER}/analysis/nw/nw

DATA=$1
#DATA=seed_10_55
#DATA=seed_56_95

rm -f ${DATA}.${PROG};
for FAMSEQ1 in `ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/${DATA}/*.fa`;
do
  for FAMSEQ2 in `ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/${DATA}/*.fa`;
  do
    echo $FAMSEQ1 > tmp.famseq
    echo $FAMSEQ2 >> tmp.famseq
    if [[ `sort tmp.famseq | head -1` == $FAMSEQ1 ]];
    then

      echo "echo $FAMSEQ1 $FAMSEQ2 >> ${DATA}.${PROG};"
      echo "{ time ${EXE} $FAMSEQ1 $FAMSEQ2; } &>> ${DATA}.${PROG};"

    fi;
  done;
done;


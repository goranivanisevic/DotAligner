BIN=/home/stesee/bin
HDOTALIGNER=/home/stesee/DotAligner
PROG=pmcomp
EXE=$BIN/pmcomp.pl

DATA=$1
#DATA=seed_10_55.new
#DATA=seed_56_95.new

rm -f ${DATA}.${PROG};
for FAMSEQ1 in `ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/${DATA}/*_dp.ps`;
do
  for FAMSEQ2 in `ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/${DATA}/*_dp.ps`;
  do
    echo $FAMSEQ1 > tmp.famseq
    echo $FAMSEQ2 >> tmp.famseq
    if [[ `sort tmp.famseq | head -1` == $FAMSEQ1 ]];
    then

      echo "echo $FAMSEQ1 $FAMSEQ2 >> ${DATA}.${PROG};"
      TMP=`sed -n 92p $FAMSEQ1 | wc -c`
      LEN1=$((TMP-2))
      TMP=`sed -n 92p $FAMSEQ2 | wc -c`
      LEN2=$((TMP-2))
      if [[ $((LEN2-LEN1)) -gt 10 || $((LEN1-LEN2)) -gt 10 ]]; then
	echo "{ time ${EXE} $FAMSEQ1 $FAMSEQ2 -D 20 | tail -5; } &>> ${DATA}.${PROG};"
      elif [[ $((LEN2-LEN1)) -gt 5 || $((LEN1-LEN2)) -gt 5 ]]; then
        echo "{ time ${EXE} $FAMSEQ1 $FAMSEQ2 -D 10 | tail -5; } &>> ${DATA}.${PROG};"
      else
        echo "{ time ${EXE} $FAMSEQ1 $FAMSEQ2 | tail -5; } &>> ${DATA}.${PROG};"
      fi

    fi;
  done;
done;


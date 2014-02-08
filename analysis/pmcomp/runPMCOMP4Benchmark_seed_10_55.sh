BIN=/home/stesee/bin
HDOTALIGNER=/home/stesee/DotAligner

rm -f seed_10_55.pmcomp;
for FAMSEQ1 in `ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/seed_10_55/*_dp.ps`;
do
  for FAMSEQ2 in `ls ${HDOTALIGNER}/analysis/rfam/selected_PIDs/seed_10_55/*_dp.ps`;
  do
    if [[ `ls -l $FAMSEQ1 | awk '{print $5}'` -gt `ls -l $FAMSEQ2 | awk '{print $5}'` ]];
    then

      echo "echo $FAMSEQ1 $FAMSEQ2 >> seed_10_55.pmcomp;"
      TMP=`sed -n 92p $FAMSEQ1 | wc -c`
      LEN1=$((TMP-2))
      TMP=`sed -n 92p $FAMSEQ2 | wc -c`
      LEN2=$((TMP-2))
      if [[ $((LEN2-LEN1)) -gt 10 || $((LEN1-LEN2)) -gt 10 ]]; then
	echo "{ time $BIN/pmcomp.pl $FAMSEQ1 $FAMSEQ2 -D 20 | tail -5; } &>> seed_10_55.pmcomp;"
      elif [[ $((LEN2-LEN1)) -gt 5 || $((LEN1-LEN2)) -gt 5 ]]; then
        echo "{ time $BIN/pmcomp.pl $FAMSEQ1 $FAMSEQ2 -D 10 | tail -5; } &>> seed_10_55.pmcomp;"
      else
        echo "{ time $BIN/pmcomp.pl $FAMSEQ1 $FAMSEQ2 | tail -5; } &>> seed_10_55.pmcomp;"
      fi

    fi;
  done;
done;


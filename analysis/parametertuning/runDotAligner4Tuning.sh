HDOTALIGNER=/home/stesee/DotAligner

for RFAM in `cat ${HDOTALIGNER}/analysis/parametertuning/rfam_pairwise_threesomes.fasta | awk 'BEGIN{count=0}/^>/{sub(">",""); if(count){print rfam":"$1;count=0}else{rfam=$1;count=1}}'`;
do
  RFAM1=${RFAM%%:*}
  RFAM2=${RFAM##*:}
  FAMSEQ1=${HDOTALIGNER}/analysis/parametertuning/data/rfam_pairwise_threesomes/${RFAM1}.pp;
  FAMSEQ2=${HDOTALIGNER}/analysis/parametertuning/data/rfam_pairwise_threesomes/${RFAM2}.pp;
  for DA_K in 0 0.1 0.3 0.5;
  do
    for DA_A in 0 -0.05 -0.1 -0.2;
    do
      for DA_B in 0 -0.05 -0.1 -0.2;
      do
        for DA_R in 0 5 10;
        do
          for DA_T in 0.25 0.5;
          do
	    for DA_L in 0 0.5;
	    do
#	      echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\tFALSE\tM\t15\tSEQALN\tTRUE' >> dotaligner.out; { time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L --pnull 0.0005 -m 15 --seqaln; } &>> dotaligner.out"
#	      echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\tFALSE\tM\t35\tSEQALN\tTRUE' >> dotaligner.out; { time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L --pnull 0.0005 -m 35 --seqaln; } &>> dotaligner.out"
	      echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\tFALSE\tM\t25\tSEQALN\tTRUE' >> dotaligner.out; { time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L --pnull 0.0005 --seqaln; }  &>> dotaligner.out"

	      for DA_S in 0 15;
	      do  
#	        echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\t'$DA_S'\tM\t15\tSEQALN\tFALSE' >> dotaligner.out; { time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L -s $DA_S --pnull 0.0005 -m 15; } &>> dotaligner.out"
#      	        echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\t'$DA_S'\tM\t35\tSEQALN\tFALSE' >> dotaligner.out; { time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L -s $DA_S --pnull 0.0005 -m 35; } &>> dotaligner.out"
	        echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\t'$DA_S'\tM\t25\tSEQALN\tFALSE' >> dotaligner.out; { time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L -s $DA_S --pnull 0.0005; } &>> dotaligner.out"

	      done;
	    done;
	  done;
        done;
      done;
    done;
  done;
done


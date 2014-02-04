HDOTALIGNER=/home/stesee/DotAligner

for RFAM in `cat ${HDOTALIGNER}/analysis/rfam/selected_PIDs/seed_56_95.fasta | awk '/^>/{split($1,a,"_"); sub(">","",a[1]); print a[1]}' | sort -u`;
do
  for FAMSEQ1 in `ls ${HDOTALIGNER}/analysis/parametertuning/data/seed_56_95/${RFAM}*.pp`;
  do
    for FAMSEQ2 in `ls ${HDOTALIGNER}/analysis/parametertuning/data/seed_56_95/${RFAM}*.pp`;
    do
      if [[ `ls -l $FAMSEQ1 | awk '{print $5}'` -gt `ls -l $FAMSEQ2 | awk '{print $5}'` ]];
      then
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
		    echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\tFALSE\tM\t15\tSEQALN\tTRUE' >> dotaligner.out"
		    echo "{ time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L --pnull 0.0005 -m 15 --seqaln; } &>> dotaligner.out"
		    echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\tFALSE\tM\t35\tSEQALN\tTRUE' >> dotaligner.out"
		    echo "{ time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L --pnull 0.0005 -m 35 --seqaln; } &>> dotaligner.out"
	            echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\tFALSE\tM\tDEFAULT\tSEQALN\tTRUE' >> dotaligner.out"
		    echo "{ time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L --pnull 0.0005 --seqaln; } &>> dotaligner.out"

		    for DA_S in 0 15;
		    do  
		      echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\t'$DA_S'\tM\t15\tSEQALN\tFALSE' >> dotaligner.out"
		      echo "{ time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L -s $DA_S --pnull 0.0005 -m 15; } &>> dotaligner.out"
      		      echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\t'$DA_S'\tM\t35\tSEQALN\tFALSE' >> dotaligner.out"
		      echo "{ time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L -s $DA_S --pnull 0.0005 -m 35; } &>> dotaligner.out"
		      echo "echo -e 'PARAMETER\t'$FAMSEQ1'\t'$FAMSEQ2'\tK\t'$DA_K'\tA\t'$DA_A'\tB\t'$DA_B'\tR\t'$DA_R'\tT\t'$DA_T'\tL\t'$DA_L'\tS\t'$DA_S'\tM\tDEFAULT\tSEQALN\tFALSE' >> dotaligner.out"
		      echo "{ time ${HDOTALIGNER}/bin/DotAligner -d $FAMSEQ1 -d $FAMSEQ2 -k $DA_K -a $DA_A -b $DA_B -r $DA_R -t $DA_T -l $DA_L -s $DA_S --pnull 0.0005; } &>> dotaligner.out"

		    done;
	          done;
	        done;
	      done;
	    done;
	  done;
	done;
      fi;
    done;
  done;
done


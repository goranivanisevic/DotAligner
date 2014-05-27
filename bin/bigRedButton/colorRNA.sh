#!/bin/bash
# Copyright (C) 2014 Stefan E Seemann <seemann@rth.dk>
if [ $# == 0 ]; then
	echo "colorRNA.sh \"((((.((....((((...(((((......).))))..))))..)..).))))\" test.fa test.pp output_prefix [hg18]"
	exit
fi

BIN_DIR=/home/stesee/bin

dt=`date --rfc-3339=ns | sed 's/[ +]/_/g'`
mkdir .f$dt

$BIN_DIR/fa2aln.pl $2 > .f$dt/test.aln
if [ $# == 5 ]
then
	CONS=`cat .f$dt/test.aln | grep $5 | awk '{print $2}'`
	TARGET=$5
else
	CONS=`$BIN_DIR/consensusSeq.pl .f$dt/test.aln`
	TARGET="hg18"
fi

head -3 .f$dt/test.aln > .f$dt/test2.aln
cat .f$dt/test.aln | grep $TARGET >> .f$dt/test2.aln

for i in `cat $2 | $BIN_DIR/seqidentity.pl - 0 | grep $TARGET | $BIN_DIR/sortColNum.pl - 3 | tac | awk -v target="$TARGET" '{if (match($1,target)) {print $2}; if (match($2,target)) {print $1}}' | xargs`; do k=${i%%:*}; cat .f$dt/test.aln | grep $k >> .f$dt/test2.aln; done
mv .f$dt/test2.aln .f$dt/test.aln

PETFOLD_ss=`cat $1`
$BIN_DIR/mycoloraln.pl -r -ss $1 -pp $3 .f$dt/test.aln > $4_aln.ps
$BIN_DIR/plotRNAstruct $CONS $PETFOLD_ss .f$dt/testss.ps
$BIN_DIR/mycolorrna.pl .f$dt/testss.ps .f$dt/test.aln > $4_ss.ps
rm -rf .f$dt RNAplot.ps

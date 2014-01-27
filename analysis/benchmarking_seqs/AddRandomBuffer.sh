#!/bin/bash

MAXLEN=250

#get sequence length
grep ">" -A 1 ${1} | grep -v -e "--" | while read line
do 
	if [[ `echo $line | cut -c 1` = ">" ]]
	then
		echo $line"_buffered_to_"$MAXLEN
	else
		LENGTH=$(echo ${line} | wc -m)
		BUFFER=$(( ${MAXLEN} -  ${LENGTH} )) 
	#	echo $BUFFER
		# echo "Seq len: "$LENGTH", Buffer to add: "$BUFFER

		#this is where to split the buffer around sequence, e.g.:
		#given SPLIT = 5 and BUFFER = 8
		#  xxxxxSSSSSSSSSSxxx  | x= random base; S = RFAM sequence
		#  12345----------678
		SPLIT=${RANDOM}
		let "SPLIT %= ${BUFFER}"
		BUFFERED=`grep ">" -A 1 ${1} | grep -v -e "--" | shuffle -d - | grep -v ">" | while read shuffled
	 	do 
	 		echo -n $shuffled
	 	done` 
		
		echo -n `echo $BUFFERED | cut -c 1-$SPLIT`
		echo -n $line
		let "SPLIT += 1"
	 	echo $BUFFERED | cut -c $SPLIT-$BUFFER
	fi
done


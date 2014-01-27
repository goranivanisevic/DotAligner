#!/bin/bash

MAXLEN=250

#get sequence length
grep ">" -A 1 ${1} | grep -v -e "--" | while read line
do 
	if [[ $line | cut -c 1 -eq ">" ]]
	then
		echo $line_buffered_to_$MAXLEN
	else
		LENGTH=$(echo ${line} | wc -m)
		BUFFER=$(( ${MAXLEN} -  ${LENGTH} )) 
		# echo "Seq len: "$LENGTH", Buffer to add: "$BUFFER

		#this is where to split the buffer around sequence, e.g.:
		#given SPLIT = 5 and BUFFER = 8
		#  xxxxxSSSSSSSSSSxxx  | x= random base; S = RFAM sequence
		#  12345----------678
		SPLIT=$(( ${RANDOM} %= ${BUFFER} ))
		BUFFERED=`grep ">" -A 1 ${1} | grep -v -e "--" | shuffle -d - | grep -v ">" | while read line
	 	do 
	 		echo -n $line 
	 	done` 
		
		echo -n $BUFFERED | cut -c 1-$SPLIT
		esho -n $1
	 	echo $BUFFERED | cut -c $SPLIT-$BUFFER
	fi
done


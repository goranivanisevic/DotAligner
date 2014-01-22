#!/bin/bash

#FILE_LIST='/home/marsmi/data/clustering/structures.list'
FILE_LIST=$1
STRUCTURES=$( wc -l ${FILE_LIST} | awk '{print $1}')
i=1
while [ $i -lt $STRUCTURES ]; do 
 	for (( j = $i+1 ; j <= $STRUCTURES; j++ )); do 
     	echo -n `head -n $i $FILE_LIST | tail -n 1`" "
     	echo -n `head -n $j $FILE_LIST | tail -n 1`" "
     	echo $i" "$j
	done
	i=$(($i+1))
done


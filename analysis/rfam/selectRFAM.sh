#!/bin/bash
for file in $* 
do 
 SIZE=$(head -n 2 $file | tail -n 1 | sed 's/\.//g' | wc -c )
 DEPTH=$(wc -l ${file} | awk '{print $1}')
 DEF=$(grep "GF DE" ${file%.fasta}| head -n 1 | cut -c 11-100 | sed 's/ /_/g')
 STRUCTURE=$(grep "GF SS" ${file%.fasta}| head -n 1 | cut -c 11-100 | sed 's/ /_/g')
 if [[ $DEPTH -gt 180 && $SIZE -gt 100 && $SIZE -lt 200 ]] 
 then 
#  cp $file ../sampled/${file%%.*}.fasta
   echo "#"${file%%.*}" "$SIZE" "$DEPTH" "$DEF" "$STRUCTURE
# GF DE/SS
 else
   echo ${file%%.*}" "$SIZE" "$DEPTH" "$DEF" "$STRUCTURE
 fi 
done

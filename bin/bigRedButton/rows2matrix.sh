#!/bin/bash
# takes a list of matrix coordinates and prints out a matrix for R
# tab delimited sorted (coord_x,coord_y)n input:
# path/to/file1 path/to/file2 coord_x coord_y score
num_cols=$(cut -d " " -f 3 $1 | grep -w 1 | wc -l)
num_cols=$(($num_cols + 1))

echo fill in the juicy bits
while read line; do 
	x=$(echo $line | cut -d " " -f 3)
	y=$(echo $line | cut -d " " -f 4)
	echo "(x,y)="$x","$y
	eval matrix_`echo $x`_`echo $y`=$(echo $line | cut -d " " -f 5)
done < $1

echo fill in the rest
i=1
j=1
while [[ $j -le $num_cols ]]; do
	while [[ $i -le $j ]]; do
		if [[ $i -eq $j ]]; then 
			eval matrix_`echo $i`_`echo $j`=1000000
		else
			eval matrix_`echo $i`_`echo $j`=$matrix_$j_$i
		fi	
		i=$(($i+1))
	done
	j=$(($j+1))
done
echo print it out
i=1
j=1
while [[ $j -le $num_cols ]]; do
	while [[ $i -lt $num_cols ]]; do
		echo -n $matrix_$i_$j
		echo -n "	"
		i=$(($i+1))
	done
	echo $matrix_$i_$j
	j=$(($j+1))
done



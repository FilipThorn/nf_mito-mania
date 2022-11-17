#!/bin/bash -l

depth=$1
name=$2

i=$(awk '{sum+=$3} END { print (sum/NR)*3}' $depth)

cut -f 3 $depth | sort| uniq | while read X; do awk -v X=$X '($3==X) { printf("%s\t%d\t%d\n",$1,$2,int($2)+1);}' $depth | sort -k1,1 -k2,2n | bedtools merge -i - | sed "s/\$/\t${X}/" ; done > mask.txt

awk '($4 < 20 || $4 > $i ) {print $1,$2,$3}' mask.txt | sort -k 1,1 -k2,2n -  > ${name}_mask_maxDP_${i}.txt

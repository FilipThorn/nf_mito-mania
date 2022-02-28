#!/bin/bash -l

#Script provided by Mario Ernst THANK YOU

echo "starting to convert coverage file $1 to bed"
tail -n+1 $1 | sort -u -nk1,1 -nk2,2 | perl -ane 'END {print " $F[1]"} next if $p[0] == $F[0] && $F[1] == $p[1] + 1; print " $p[1]\n@F"; } continue { @p = @F;' > $1.bed
​
sed -i -e 's/ /\t/g' $1.bed

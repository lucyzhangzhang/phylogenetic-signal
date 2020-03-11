#!/bin/bash

#converts phylip alignment files to nexus files

OUT=$1.nexus

if [[ -e $OUT ]]; then
    rm $OUT
fi


echo "#nexus" >> $OUT

echo "begin data;" >> $OUT

tax=$(head -n 1 $1 | awk '{print $1}')
char=$(head -n 1 $1 | awk '{print $2}')

echo -e "\tdimensions ntax=$tax nchar=$char;" >> $OUT

#does account for RNA
if [[ $(sed '1d' $1 | cut -d" " -f2- | grep -ie '[BDEFHIJKLMNOPQRSVW]') ]]; then
    type="protein"
else 
    type="dna"
fi

echo -e "\tformat datatype=$type missing=? gap=-;" >> $OUT

echo -e "\tmatrix" >> $OUT

lineNum=$(sed '1d' $1 | wc -l)
narg=$(echo "scale=0; $lineNum/$tax"| bc)

for ((i=1;i<=$tax;i++)) ; do
    echo -ne "\t\t" >> $OUT
#    printf '%s' $(sed '1d' $1 | sed '/^$/d' | sed -n $i'~'$tax'p' | perl -pe 's/\s/\t/') >> $OUT
    out=$(printf '%s' $(sed '1d' $1 | sed -n $i'~'$tax'p'))
    echo $out | perl -pe 's/(.{10})/$1\t/' >> $OUT
done

echo -e "\t;" >> $OUT 

echo "end;" >> $OUT

#!/bin/bash

#lmao
# usage: script dir E.fa A.fa

cd $1

if [ -d tmp ]
then
    rm -rf tmp
    mkdir tmp
fi

cd tmp 

grep '>' $2 > E

sed -i 's/>//g' E

while read line
do
done < E

rm -f E 

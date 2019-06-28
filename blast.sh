#!/bin/bash
# Usage: blast <transcript folders>
DIR=$PWD

if [[ ! -d $DIR/blast ]]; then
    mkdir $DIR/blast
fi

#Recognize .fasta, .fna, and .fa files
function blast_hk () {
    if [[ -e $DIR/blast/$3.out ]]; then
        rm $DIR/blast/$3.out
    fi
    for i in `ls $2/housekeep/*.f*a*`; do
    blastn -query $i -db $1/*.f*a -outfmt 6 -evalue 0.00001 -max_hsps 1 -max_target_seqs 10 -task "dc-megablast" | sort -r -k4,4n  >> $DIR/blast/$3.out
    done
    for i in `ls $2/housekeep/transcripts/*.f*a*`; do
    blastn -query $i -db $1/*.f*a -outfmt 6 -evalue 0.00001 -max_hsps 1 -max_target_seqs 10 -task "dc-megablast" >> $DIR/blast/$3.out
    done
}

#Make blast database
function makedb () {
    if ls $1 | grep '(nhr|nin|nsq)'; then
        echo "Database exists"
    else    
        makeblastdb -dbtype nucl -in $1/*.f*a -out $1/*.f*a
    fi
}

while read spec; do
    DATA=$DIR/transcripts/$spec
    makedb $DATA
    blast_hk $DATA $DIR/genes $spec
done < $1

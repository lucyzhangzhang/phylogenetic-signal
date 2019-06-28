#!/bin/bash
# Usage: blast <transcript folders>
DIR=$(dirname "$1")

if [[ ! -d $DIR/blast ]]; then
    mkdir $DIR/blast
fi

#Recognize .fasta, .fna, and .fa files
function blast_genes () {
    for i in $(ls $2/housekeep/*.f*a* & ls  $2/phos/*.f*a); do
    blastn -query $i -db $1/*.f*a -outfmt 6 -evalue 0.00001 -max_hsps 1 -max_target_seqs 1 -task "dc-megablast" | sort -nrk 4,4  >> $DIR/blast/dcMega.out
    blastn -query $i -db $1/*.f*a -outfmt 6 -evalue 0.00001 -max_hsps 1 -max_target_seqs 1 | sort -nrk 4,4  >> $DIR/blast/blast.out
    done
}

#Make blast database
function makedb () {
    if $(ls $1 | grep '(nhr|nin|nsq)'); then
        echo "Database exists"
    else    
        makeblastdb -dbtype nucl -in $1/*.f*a -out $1/*.f*a
    fi
}

while read spec; do
    echo -e "=============\n====$spec====\n=============" | tee -a $DIR/blast/dcMega.out $DIR/blast/blast.out
    DATA=$DIR/transcripts/$spec
    makedb $DATA
    blast_genes $DATA $DIR/genes
done < $1

#!/bin/bash
# Usage: transcriptSearch.sh <configfile>

DIR=$(dirname "$1")

function getTranscript () {
    #Gene ID, gene locus and organism name
    ID=$1
    LOC=$2
    ORG=$3

    #Handy sed! | remove last line, as it matches the beginning of the next fasta
    sed -n "/$LOC/,/>/p" /home/lucy/scratch/physig/transcripts/$ORG/*.f*a | sed '$d' > /home/lucy/scratch/physig/genes/$ORG.$ID.$LOC.fa
}

export -f getTranscript

for spec in `cat $1`; do
    echo "getting $spec transcripts"
    parallel -j 10 -C '\t' --group getTranscript :::: $DIR/genes/$spec.id ::: $spec
    echo "get $spec finised"
done

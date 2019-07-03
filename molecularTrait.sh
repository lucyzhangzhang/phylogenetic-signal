#!/bin/bash

# Usage moleculatTrait.sh $PATH/<configFile>
# Uses the $spec.id files pasted from transcriptSearch.sh
# Parses gff files from Ensemble

DIR=$(dirname "$1")

function exonNum () {
    echo "Exon number"
    #count of all the exon entries
    # $1=IDfile, $2=GFF, $3.spec
    
    while read gene ID; do
        if [[ $3 == "Zmays" ]]; then
            grep "transcript:$ID" $2 | grep -c 'exon' >> $3.exons
        elif [[ $3 == "Esal" ]]; then
            TrID=$(grep "$ID" $DIR/transcripts/Esal/*.gtf | grep 'exon' | perl -ne '/.*transcript_id "(.*?)".*/ && print("$1\n")' | head -n 1)
            grep "$ID" $DIR/transcripts/Esal/*.gtf | grep 'exon' | grep -c "$TrID" >> $3.exons
        elif [[ $ID == "MtPIDL1" ]] && [[ $3 == "Mtrunc" ]]; then
            echo "1" >> $3.exons
        else
            grep "$ID" $2 | grep -c 'exon' >> $3.exons
        fi
    done < $1
}

function ORFLen () {
    echo "ORF length"
    #sum of all the lengths of the CDS

    while read gene ID; do
        if [[ $ID == "MtPDIL1" ]] && [[ $3 == "Mtrunc" ]]; then
            echo ""
        else
            grep $ID $2 > tmp
            if [[ $(awk '$3=="CDS"' tmp) ]]; then
                awk '$3=="CDS"' tmp | awk '{t=$4-$5;print t;}' | sed 's/-//g' | awk '{n+=$1} END {print n}' >> $3.ORF
            else
                awk '$3=="exon"' tmp | awk '{t=$4-$5;print t;}' | sed 's/-//g' | awk '{n+=$1} END {print n}' >> $3.ORF
            fi
        fi
    done < $1
    rm tmp
}

function mRNAlen () {
    echo "transcript length"
    #length of the mRNA transcript

    while read gene ID; do
        FILE=$DIR/alignment/$gene.fa

        sed -n "/$ID/,/>/p" $FILE | sed '$d' | wc -c >> $3.mRNA
    done < $1
}

function GC () {
    echo "GC content"
    #percentage of G and C for transcript

    while read gene ID; do
        FILE=$DIR/alignment/$gene.fa

        sed -n "/$ID/,/>/p" $FILE | sed '$d' | sed '/>/d' > seq
        AT=$(grep -o '[AT]' seq | wc -l) 
        GC=$(grep -o '[GC]' seq | wc -l)

        echo -e "scale=5; $GC/($GC+$AT)" | bc >> $2.GC
        
        rm seq
    done < $1                                      
}

if [[ -e molecularTrait.out ]]; then
    rm molecularTrait.out
fi

OUT=molecularTrait.out

echo -e "Org\tGene\texonNum\tORFLen\tmRNAlen\tGC" >> $OUT

sed '/^$/d' $DIR/blast/select | awk '{if($0 ~ /^#.*/) name=$0} {sub(/#/,"",name)}  /^#/{close(file);file=$NF} /./{print > name".id"}'

while read spec; do
    echo "Getting traits for $spec"

    GFF=$DIR/transcripts/$spec/*.gff*

    sed -i '/#/d' $spec.id
    IDfile=$spec.id

    exonNum $IDfile $GFF $spec &
    ORFLen $IDfile $GFF $spec &
    mRNAlen $IDfile $GFF $spec &
    GC $IDfile $spec &

    start=1
    end=$(wc -l < $spec.id)

    for ((i=start; i<=end; i++)); do
        echo "$spec" >> $spec.name
    done
    
    wait 

    paste $spec.name $IDfile $spec.exons $spec.ORF $spec.mRNA $spec.GC > $spec.trait
done < $1

cat *.trait >> $OUT

echo "Removing tmp files"
rm *.trait *.exons *.ORF *.mRNA *.GC *.id *.name

echo "Completed"

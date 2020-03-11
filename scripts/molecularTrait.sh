#!/bin/bash

# Usage moleculatTrait.sh $PATH/<configFile>
# Uses the $spec.id files pasted from transcriptSearch.sh
# Parses gff files from Ensemble

DIR=$(dirname "$1")

function getTranscript () {
    echo "Getting transcripts..."

    #get transcripts of each gene
    FA=/home/lucy/scratch/physig/transcripts/$2/*.f*a
    while read gene ID; do
        if [[ $ID == "MtPDIL1" ]] && [[ $2 == "Mtrunc" ]]; then
            cat genes/Mtrunc.IPS1.fa >> $gene.fa
        elif [[ $2 == "Zmays" ]] && [[ $gene == "IPS1" ]]; then
            cat genes/Zmays.Zm00001d022669_T001.fa >> $gene.fa
        elif [[ $2 == "Esal" ]] && [[ $gene == "IPS1" ]]; then 
            cat genes/Esal.IPS1.XM_024159506.1.fa >> $gene.fa
        else
            sed -n "/$ID/,/>/p" $FA | sed '$d' >> $gene.fa
        fi
   done < $1

}

function exonNum () {
    echo "Exon number"
    #count of all the exon entries
    # $1=IDfile, $2=GFF, $3.spec
    
    while read gene ID; do
        if [[ $3 == "Zmays" ]]; then
            grep "transcript:$ID" $2 | grep -c 'exon' >> $3.exons
#        elif [[ $3 == "Esal" ]]; then
#            TrID=$(grep "$ID" $DIR/transcripts/Esal/*.gtf | grep 'exon' | perl -ne '/.*transcript_id "(.*?)".*/ && print("$1\n")' | head -n 1)
#            grep "$ID" $DIR/transcripts/Esal/*.gtf | grep 'exon' | grep -c "$TrID" >> $3.exons
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
#        if [[ $gene == "IPS1" ]]; then
#            echo "dunno" >> $3.ORF
#        elif [[ $3 == "Esal" ]]; then
            #Eutrema is fking special asdfasdfjasdf
#            TrID=$(grep "$ID" $DIR/transcripts/Esal/*.gtf | grep 'exon' | perl -ne '/.*transcript_id "(.*?)".*/ && print("$1\n")' | head -n 1)
#            grep "$TrID" $DIR/transcripts/Esal/*.gtf | grep 'exon' | awk '{t=$4=+$5;print t}' | sed 's/-//g' | awk '{n+=$1} END {print n}' >> $3.ORF
#        else
        grep $ID $2 > tmp
            if [[ $(awk '$3=="CDS"' tmp) ]]; then
                awk '$3=="CDS"' tmp | awk '{t=$4-$5;print t;}' | sed 's/-//g' | awk '{n+=$1} END {print n}' >> $3.ORF
            else
                awk '$3=="exon"' tmp | awk '{t=$4-$5;print t;}' | sed 's/-//g' | awk '{n+=$1} END {print n}' >> $3.ORF
            fi
#        fi
    done < $1

    rm tmp
}

function mRNAlen () {
    echo "transcript length"
    #length of the mRNA transcript

    while read gene ID; do
        sed -n "/$ID/,/>/p" $gene.fa | sed '$d' | wc -c >> $3.mRNA
    done < $1
}

function GC () {
    echo "GC content"
    #percentage of G and C for transcript

    while read gene ID; do
        FILE=$gene.fa

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

#remove existing sequence files
rm *.fa 

OUT=molecularTrait.out

echo -e "Org\tGene\tID\texonNum\tORFLen\tmRNAlen\tGC" >> $OUT

sed '/^$/d' $DIR/blast/select | awk '{if($0 ~ /^#.*/) name=$0} {sub(/#/,"",name)}  /^#/{close(file);file=$NF} /./{print > name".id"}'

while read spec; do
    echo "Getting traits for $spec"

    GFF=$DIR/transcripts/$spec/*.gff*

    sed -i '/#/d' $spec.id
    IDfile=$spec.id
    
    getTranscript $IDfile $spec

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

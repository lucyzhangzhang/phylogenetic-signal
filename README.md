# Phylogenetic Signal
## Co-evolving genes AKA genes related to lncRNA and IPS1 Transcripts
- PHO2
- PHR1
- WRKY75
- SIZ1
- PHT1
- RNS1
- SULTR1
- PHF

## Housekeeping genes
- EF1a
- Cyclophilin
- Hsp90-1

## Organisms
- E. salsugineum
- A. thaliana
- O. sativa
- B. rapa
- C. rubella
- B. napus
- M. trunculata
- S. lycopersicum
- Z. Mays (Outgroup)

## Steps:
1. Download transcriptomes from SRA
2. Download genomes from NCBI or JGI
3. Align transcripts to genomes
    a) Get annotation files from ensemble
    b) Get transcript files from ensemble
4. Create annotation with existing annotation (RSEM?)
5. Extract transcript sequences

    a) Blast transcript sequences known in E. salsugineum with other transcripts

    ```ls -1 transcripts > config && blast.sh $PWD/config```

    b) Parse top hit transcript
    c) Manual go through each hit D:
        i) alternatively, somehow get ThaleMine's Phytozone homologs working (API?)
6. Align transcript sequences

    ```mafft --phylipout --nuc input > output```

7. Generate tree from alignment - Using RaxML

    ```raxmlHPC-PTHREADS -s IPS1.phy -n rax1 -m GTRCAT -o Pp3c11_223 -N 1000 -x 51 -p 51 -k -d -T 5```

8. Determine molecular traits of each group: ORF length, GC content, # of exons, transcript length
9. Calculate phylogenetic signal from R package "phylosignal"

## Notes:
* how do we know that these housekeeping genes don't interact with lncRNAs?

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
    - 18s rRNA
    - Cyclophilin
    - Hsp20.2

## Organisms:
    - E. salsugineum
    - A. thaliana
    - O. sativa
    - B. rapa
    - B. napus
    - M. trunculata
    - C. reinhardti
    - S. lycopersicum
    - P. patens

## Steps:
    1. Download transcriptomes from SRA
    2. Download genomes from NCBI or JGI
    3. Align transcripts to genomes
    4. Create annotation with existing annotation (RSEM?)
    5. Extract transcript sequences
        a) Blast transcript sequences known in E. salsugineum with other transcripts
        b) Parse top hit transcript
    6. Align transcript sequences
        mafft --phylipout --nuc input > output
    7. Generate tree from alignment - Using RaxML
        raxmlHPC-PTHREADS -s IPS1.phy -n rax1 -m GTRCAT -o Pp3c11_223 -N 1000 -x 51 -p 51 -k -d -T 5
    8. Calculate phylogenetic signal from R package "phylosignal"

## Notes:
*how do we know that these housekeeping genes don't interact with lncRNAs?

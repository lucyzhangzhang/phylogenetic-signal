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

    ```#Parsinomious tree
        ls *.fa.phy | xargs -i raxmlHPC -f d -m GTRCAT -s {} -n {}.out -p 51
        #Bootstrapping*
        ls *.phy | xargs -P 6 -i raxmlHPC -f d -m GTRCAT -s {} -N 1000 -b 51 -p 51 -n {}.bt*
        #Consensus tree
        ls *.phy | xargs -P 6 -i raxmlHPC -f b -m GTRCAT -s {} -n {}.cons -z *.boopstrap.{}.bt
    ```

8. Determine molecular traits of each group: ORF length, GC content, # of exons, transcript length
9. Calculate phylogenetic signal from R package "phylosignal"

## With regards to *M. truncatula*
* Conserved region for IPS2 found but not IPS1

    a) Possible incomplete genome assembly?

## With regards to *S. lycopersicum*
* Very short transcript (180 bp) mapped to region with 22 bp conserved region, instead usual ~520 bp

    a) Possible incomplete transcriptome assembly?

## Notes:
* how do we know that these housekeeping genes don't interact with lncRNAs?

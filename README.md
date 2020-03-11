# Phylogenetic Signal
## Co-evolving genes AKA genes related to lncRNA and IPS1 Transcripts
- PHO2
- PHR1
- WRKY75
- SIZ1
- PHT1
- RNS1
- SULTR2;1
- PHF
- PHOX (*C. reinhardtii*, putative phosphatase)

## Housekeeping genes
- EF1a
- EF1b
- EF4a
- Cyclophilin
- Hsp90-1
- Actin
- Tubulin
- Glyceraldehyde-3-phosphate dehydrogenase
- SuCoA (Succinol Co-enzyme A)
- FtSH protease
- 18S (both plants and algae)
- GADPH
- 25S rRNA (Plants only)
- 28S rRNA (*C. reinhardtii*)

To add:
- actually other housekeeping genes

## Sulfur metabolism gemes
- SULTR
- SNRK
- ???SAC1 (*C. reinhardtii*)
- ???SLT (*C. reinhardtii*)
- ???LPB1 (*C. reinhardtii*)
- SULP
- ATS
- SBP
- APK 

## Organisms
- Eutrema salsugineum
- Arabidopsis thaliana
- Oryza sativa
- Brassica rapa
- Capsella rubella
- Brassica napus
- Medicago trunculata
- Solanum lycopersicum
- Zea Mays (Outgroup)

To add:
- Amborella trichopoda
- Chlamydomonas reinhardtii (new outgroup)
- Solanum tuberosum
- Selaginella moellendorffii

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

    c) Manual go through each hit D:.i) alternatively, somehow get ThaleMine's Phytozone homologs working (API?)

6. Align transcript sequences

    ```mafft --phylipout --nuc input > output```

7. Generate tree from alignment using conserved genes - Using RaxML

    ```#Parsinomious tree
        ls *.fa.phy | xargs -i raxmlHPC -f d -m GTRCAT -s {} -n {}.out -p 51
        #Bootstrapping*
        ls *.phy | xargs -P 6 -i raxmlHPC -f d -m GTRCAT -s {} -N 1000 -b 51 -p 51 -n {}.bt*
        #Consensus tree
        ls *.phy | xargs -P 6 -i raxmlHPC -f b -m GTRCAT -s {} -n {}.cons -z *.boopstrap.{}.bt
    ```

8. Infer species trees from bootstrap trees generated from conserved genes (consensus tree is a good way to go)
    
    a) Don't forget to test for nuclotide model using `jModelTest`

9. Determine molecular traits of each group: ORF length, GC content, # of exons, transcript length
10. Calculate phylogenetic signal from R package "phylosignal"

## With regards to *M. truncatula*
* Conserved region for IPS2 found but not IPS1

    a) Possible incomplete genome assembly?

## With regards to *S. lycopersicum*
* Very short transcript (180 bp) mapped to region with 22 bp conserved region, instead usual ~520 bp

    a) Possible incomplete transcriptome assembly?

## Notes:
* how do we know that these housekeeping genes don't interact with lncRNAs?
* Are there families of genes close together that modulate Pi homeostasis such as the PHO regulon in bacteria? (The answer is likely: no)
* \*Update\* Now that I can use Phytozome API, I can find homologs (already annotated)
* Genes that don't have homologs across all the species then what?
* At least I have to build a topology from conserved genes among the 

## Phosphate starvation response *C. reinhadrtii*
- PSR1 (conserved MYB transcription factor)

## Part 2: API Scripting
How to use the Python API script

```
phytozome.py -h
```

This is a version made for phylogenetic signal analysis

Functions to implement:
* Finding all the homologous genes (as much as possible)
* Find the transcript sequences (required for alignment)
* Possible? Number of exons
* Counting GC%
* ORF length (probably through counting the CDS
* Transcript length (relatively EZ)
* Will have to keep track of ones that don't have any homologs...(which is a lot)
* GO analysis to determine what role the transcript has

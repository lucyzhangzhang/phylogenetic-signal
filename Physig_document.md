# Phylogenetic Signal
The goal of phylogenetic signal is to see whether certain traits of genes are conserved across different species or they change based on how related each species is to each other.

For this project, we want to see what phylogenetic signal the following traits have on housekeeping, phosphate-starvation-related, and all other genes:
- GC%
- Transcript length
- CDS (coding sequence) length
- Exon number

The GC% is easy to find as it's just a matter of counting the number of G's and C's in a gene and dividing by the total gene length.

The Transcript length is also easy because it's just a matter of counting how many characters are in the transcript length.

The CDS length and exon number are much more difficult as it requires the annotation file (which contains the exon and CDS information).

The phylogenetic signal will be performed in the following species (all of these species can be found in the Phytozome database):
- Eutrema salsugineum
- Arabidopsis thaliana
- Oryza sativa
- Brassica rapa
- Capsella rubella
- Brassica napus
- Medicago trunculata
- Solanum lycopersicum
- Zea Mays
- Amborella trichopoda
- Solanum tuberosum
- Selaginella moellendorffii
- Chlamydomonas reinhardtii (outgroup)

## Phylogenetic signal steps
The focus of the analysis is on *Eutrema salsugineum* so we want to find homologs of all the genes in our speices and compare them to other species.
1. Find homologs of all the genes in *E. salsugineum* in the other species
2. Based on the GO terms (if exists), determine if the genes are housekeeping, phosphate-stress-related, or "other" group.
3. Find transcript sequences of all the genes+homologes
4. Keep track of genes that don't have any homolog results, (we want to run BLAST on these genes)
5. Put the transcript sequences of each gene and their homologs in their own multi FASTA file
6. We will align these genes and determine the maximum likelihood tree for each of these FASTAs by running MAFFT and RAxML

## For Tianhe
The Phytozome API script will be revised to incorporate the following new features:
- Instead of only looking in *Arabidopsis* for homologs, you will extend the search of homologs in all of the species listed above. I.e. change the argParse of "A. thaliana" to instead a list of all the organisms, perhaps the arg will be an option to read a file with all the organisms listed:
- Same thing will happen for looking for FASTA files, you will combined the FASTA of homologs from one gene into one file by appending them to teach other. I.e the FASTA for a gene called RNS1 will have:
```
> Eutrema RNS1
XXXXXXXXXXX
> Arabidopsis RNS1
XXXXXXXXXXX
> O. sativa RNS1
XXXXXXXXXXX
...
> C. reinhardtii RNS1
XXXXXXXXXXX
```
- The gene list of all *Eutrema* genes is here:
```
/home/lucy/Eutrema/FPKM/droughtName.uniq
```
- The entire transcripts database of *Eutrema* genes is here:
```
/2/scratch/lucy/misc/2transcripts.fa
```
- For the transcripts that don't have any homologs, you can extract them from the transcripts database like so using a program called `samtools`:
```
samtools faidx 2transcripts.fa <transcript_name> -o transcript_name.fa
```
This runs the FASTA indexing tool on 2transcripts.ca to extract the single FASTA of the name transcript\_name and puts it in a file called transcript\_name.fa. Try it out!!

#!/usr/bin/env python3
# need to have Intermine installed
# you can do a local install through conda

from intermine.webservice import Service  # intermine Phytomine API service
import argparse # reading command line arguments
import os # reading input from the STDIN
import re # required for printing fasta format

usage = "usage %prog GeneList \"OrganismShortName\""

parser = argparse.ArgumentParser(
        description = "Interacts with the Phytozome API, outputs two files: the matches (matches.tab) and a multifasta of the mismatches (unmatched.fa)"
        )
parser.add_argument("-l", "--list", dest = "geneList", 
        required = True, metavar = "Gene List", 
        # in the case of phylogenetic signal, that would be everything, as much as possible lol
        help = "List of names of genes you want to search for, each on their own line")
parser.add_argument("-n", "--name", dest = "Sname", 
        required = True, type = str,
        # there's a random newline here for some reason
        metavar = "Organism list", 
        help = "list of Organisms to look for homologs in")
parser.add_argument("-o", "--outDir", dest = "outDir",
        required = False, metavar = "Output directory location",
        type = str, help = "Location to put the output files")

service = Service("https://phytozome.jgi.doe.gov/phytomine/service")

args = parser.parse_args()

#path of the input file
fpath = os.path.dirname(args.geneList)

# test if outDir is given
if args.outDir is None:
    out = fpath
else:
    out = args.outDir

# remove if exists
def deletfile(filepath):
    if os.path.exists(filepath):
        try:
            os.remove(filepath)
        except:
            print("error when deleting file", filepath)

outFiles = [os.path.join(out, "matches.tab"), os.path.join(out, "unmatched.fa")]
[deletfile(i) for i in outFiles]

# appends to existing file instead of overwriting with "w"
matches=open(os.path.join(out, "matches.tab"), "a")
unmatched=open(os.path.join(out, "unmatched.fa"), "a")

with open(args.geneList) as gList:
    for gName in gList:
        gName = gName.rstrip()
        print(gName)

        # initiate new query
        query = service.new_query("Homolog")

        # specifying query parameters
        query.add_view(
            "gene.primaryIdentifier", "ortholog_gene.primaryIdentifier",
            "ortholog_organism.shortName", "relationship"
        )

        # query specifications
        query.add_constraint("gene.primaryIdentifier", "CONTAINS", str(gName), code = "A")
        query.add_constraint("organism.shortName", "=", args.Sname, code = "B")
        query.add_constraint("ortholog_organism.shortName", "=", "E. salsugineum", code = "C")

        if len(query) == 0:
            query2 = service.new_query("Gene")
            query2.add_view("primaryIdentifier", "transcripts.primaryIdentifier",
                    "transcripts.sequence.residues"
            )
            query2.add_sort_order("Gene.transcripts.primaryIdentifier", "ASC")
            query2.add_constraint("transcripts.primaryIdentifier", "CONTAINS", gName, code = "A")
            query2.add_constraint("organism.shortName", "=", args.Sname, code = "B")
            for row in query2.rows():
                faFormat=re.sub("(.{75})", "\\1\n", row["transcripts.sequence.residues"], 0, re.DOTALL)
                unmatched.write(">" + row["primaryIdentifier"] + "\n" + faFormat + "\n")
#                print(">" + row["primaryIdentifier"] + "\n" + faFormat)
        else:
            for row in query.rows():
                matches.write(row["gene.primaryIdentifier"] + "\t" + 
                    row["ortholog_gene.primaryIdentifier"] + "\t" +
                    row["ortholog_organism.shortName"] + "\t" +
                    row["relationship"] + "\n")
#                print(row["gene.primaryIdentifier"], 
#                    row["ortholog_gene.primaryIdentifier"],
#                    row["ortholog_organism.shortName"],
#                    row["relationship"])
#        for row in query.results("tsv"):
#            print(row)


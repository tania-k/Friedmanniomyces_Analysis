#!/usr/bin/env python3

import csv, argparse, os, sys, re
import Bio
from Bio.Seq import Seq
from Bio import SeqIO


# we will use argparse to support command line arguments
# see https://docs.python.org/3/library/argparse.html

parser = argparse.ArgumentParser(description='Get Duplicated Ortholog sets from OrthoFinder results.')
parser.add_argument('--input_folder','--input',
                    default = 'input/OrthoFinder/Results_Jan16',
                    help='OrthoFinder Results folder')

parser.add_argument('--db','--cds_db',
                    default = 'Dothid.cds.fa',
                    help='CDS sequence database in fasta format')

parser.add_argument('--out','--cds_out',
                    default = 'CDS_aligned_per_species',
                    help='folder to write CDS sequences split out by species')

args = parser.parse_args()
cds_seq_dict = SeqIO.index(args.db, "fasta")

if not os.path.isdir(args.out):
    os.mkdir(args.out)

# we need to remove the prefix because these are not consistent
# between the CDS file and the protein files used for searching.  In
# my best practices I have a prefix so it is easy to know which
# species a gene is from, but we can infer that already from the LOCUS
# anyways.

orthtable = os.path.join(args.input_folder,'Orthogroups','Orthogroups.tsv')
with open(orthtable,"r") as fh:
    parseOrth = csv.reader(fh,delimiter="\t")
    header = next(parseOrth)
    species_folders = []
    for h in header[1:]:
        species = re.sub(r'\.aa$','',h)
        species_folders.append(os.path.join(args.out,species))
        if not os.path.isdir(species_folders[-1]):
            os.mkdir(species_folders[-1])

    for row in parseOrth:
        og = row[0]
        col = 0
        for sp in row[1:]:
            outdir = species_folders[col] # get the output folder name
                                          # based on the order of headers
            genes = sp.split(", ")
            if len(genes) == 2: # only keep pairwise duplicates for now
                seqs = []
                pepseqs = []
                for g in genes:
                    genename = re.sub(r'^([^\|]+)\|','',g)
                    if genename in cds_seq_dict:
                        cdsseq = cds_seq_dict[genename]
                        pepseq = cdsseq.translate(to_stop=True)
                        pepseq.id = genename
                        pepseq.description = cdsseq.description
                        if len(pepseq) + 2 < (len(cdsseq) / 3):
                            print("Skipping truncated/pseudogene {}".format(genename))
                            continue
                        pepseqs.append(pepseq)
                        seqs.append(cdsseq)
                    else:
                        print("cannot find {} ({}) in database".format(genename,g))
                SeqIO.write(seqs,
                            os.path.join(outdir,"{}.cds.fas".format(og)),
                            'fasta')
                SeqIO.write(pepseqs,
                            os.path.join(outdir,"{}.pep.fas".format(og)),
                            'fasta')

            col += 1

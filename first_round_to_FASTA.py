#/bin/python
#
# Using a list of protein IDs of HMMER-predicted proteins
# and your directory of protein FASTAs, 
# pull out protein sequences for the 
# second pass with HMMER as a FASTA file
#

from Bio import SeqIO
import os, fnmatch, re, glob


# Place input files here. Change for different MAMPs
hits = []
ID_list = 'results/hmmer/ef_tu/first_hits.txt'

protein_FASTA_directory = 'proteins/*/*/*.faa'

output_fasta = open('results/hmmer/ef_tu/first_round_candidates.faa', 'w')


# Open the table output from HMMER
with open(ID_list) as hit_list:
  for line in hit_list:
    hits.append(line.strip())
    

# Get all the protein FASTA files in a directory
# and print any matching protein IDs 
# (I use glob because I organized FASTA files 
# in folders by taxonomic class)
for filename in glob.glob(protein_FASTA_directory):
  filename_without_extension, extension = os.path.splitext(filename)
  with open(filename) as protein_fasta:
    for record in SeqIO.parse(protein_fasta, "fasta"):
        for hit in hits:
          if hit == record.id:
            output_fasta.write(record.format("fasta"))
            
output_fasta.close()
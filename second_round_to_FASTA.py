#/bin/python
#
# Using a list of protein IDs of HMMER-predicted proteins
# and your directory of protein FASTAs, 
# pull out protein sequences as a FASTA file
# and count how many there are per genome
#

from Bio import SeqIO
import os, fnmatch, re, glob

# Add input files here. Change for different MAMPs 
hits = []
hits_no_dupes = []

ID_list = 'results/hmmer/ef_tu/final_hits.txt'

protein_fasta_directory = 'proteins/*/*/*.faa'

output_fasta = open("results/hmmer/ef_tu/final_candidates.faa", "w")


# Create an array of protein IDs without duplicates 
with open(ID_list) as hit_list:
  for line in hit_list:
    # Use an index of 11 for the elfs, 19 for the csps, and 11 for flgs.
    hits.append(line.split(' ')[11])
    
for hit in hits:
  if hit not in hits_no_dupes:
    hits_no_dupes.append(hit)
    
    
# Write a FASTA file of all MAMP candidates by 
# matching the protein ID to the corresponding sequence
# in each genome
for filename in glob.glob(protein_fasta_directory):
  filename_without_extension, extension = os.path.splitext(filename)
  with open(filename) as protein_fasta:
    for record in SeqIO.parse(protein_fasta, "fasta"):
        for hit in hits_no_dupes:
          if hit == record.id: 
            record.id = filename_without_extension.split('/')[3][0:15].strip()+"_"+hit
            record.description = ""
            output_fasta.write(record.format("fasta"))
            
output_fasta.close()


# From a list of protein IDs, make files of all the flg sequences per genome. 
# Make a tab-delimited file of how many flg sequences per genome. 
EFTus_per_genome = {}
EFTu_count = 0


for filename in glob.glob('proteins/*/*/*converted.faa'):
  filename_without_extension, extension = os.path.splitext(filename)
  with open(filename) as protein_fasta:
    EFTu_count = 0
    output=open(filename_without_extension+"_EFTus.faa", "w")
    for record in SeqIO.parse(protein_fasta, "fasta"):
      for hit in hits_no_dupes:
        if hit == record.id:
          output.write(record.format("fasta"))
          EFTu_count += 1
    EFTus_per_genome[filename_without_extension.split('/')[3][0:15].strip()] = filename_without_extension.split("/")[2]+"\t"+str(EFTu_count)
    output.close()
    
EFTu_count_tsv = open("results/hmmer/ef_tu/EFTu_counts.tsv", "w")
EFTu_count_tsv.write("genome\tgenus\tnumber_of_EFTus\n")
for genome in EFTus_per_genome:
  EFTu_count_tsv.write(genome+"\t"+str(EFTus_per_genome[genome])+"\n")
EFTu_count_tsv.close() 
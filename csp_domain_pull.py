#bin/python

# Given a TSV of WP tags,
# the start location of the alignment, 
# and the end location of the alignment, 
# produce a FASTA file with these domains
# and a TSV file with metadata of what bacteria these
# domains belong to.

import os, glob
from Bio import SeqIO

# Input file: tsv with WP tags and loci
# Input directory: points to annotated protein files
# in separate folders based on bacterial family 

input_file = "results/hmmer/cold_shock_prot/domain_extraction/csp_prot_domain_loci.txt"

domains = [] 

with open(input_file, "r") as f:
  for line in f:
    domains.append(line.strip())
    
domain_seqs = {}

output_fasta = "results/hmmer/cold_shock_prot/domain_extraction/extracted_csp_domains.fasta"
output_metadata = "results/hmmer/cold_shock_prot/domain_extraction/WP_metadata.tsv"
    
# Loop through list of domains
# if there's a WP tag match to one of the protein files,
# record the strain as designated in the proper directory 
# and grab the corresponding domain 

for filename in glob.glob("proteins/*/*/*protein.faa"):
  filename_without_extension, extension = os.path.splitext(filename)
  with open(filename) as protein_fasta:
    for record in SeqIO.parse(protein_fasta, "fasta"):
      for domain in domains:
        WP_tag = domain.split(" ")[0]
        start = domain.split(" ")[1]
        end = domain.split(" ")[2]
        if record.id == WP_tag:
          if record.id not in domain_seqs and [str(record.seq)[int(start)-1:int(end)-1], filename_without_extension.split("/")[2], start, end] not in domain_seqs.values():
            domain_seqs[record.id] = [str(record.seq)[int(start)-1:int(end)-1], filename_without_extension.split("/")[2], start, end]
            print([str(record.seq)[int(start)-1:int(end)-1], filename_without_extension.split("/")[2], start, end])
          elif [str(record.seq)[int(start)-1:int(end)-1], filename_without_extension.split("/")[2], start, end] not in domain_seqs.values():
            domain_seqs[record.id+"_ex"] = [str(record.seq)[int(start)-1:int(end)-1], filename_without_extension.split("/")[2], start, end]
            print([str(record.seq)[int(start)-1:int(end)-1], filename_without_extension.split("/")[2], start, end])
      
# Write outputs
# One will be a FASTA file to use to generate an alignment
# One will be a TSV containing domain metadata to use when building the csp domain tree
with open(output_fasta, "w") as fasta_output:
  for key, value in domain_seqs.items():
    fasta_output.write(">"+key+"\n")
    fasta_output.write(value[0]+"\n")
  
with open(output_metadata, "w") as metadata_output:
  metadata_output.write("WP_tag\tgenus\tdomain_seq\tstart\tend\n")
  for key, value in domain_seqs.items():
    metadata_output.write(key+"\t"+value[1]+"\t"+value[0]+"\t"+value[2]+"\t"+value[3]+"\n")
    
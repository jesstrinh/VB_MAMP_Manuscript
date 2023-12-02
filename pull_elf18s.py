#/bin/python
#
# Given a FASTA file of MAMP candidates
# and a table of GenBank accession numbers and species names
# pull out epitope (ex: elf18) sequences from your dataset 
# and create summary tables of all sequences and unique sequences
#
# NOTE: This may not capture ALL possible elf18 sequences. 
# Manually-curate any 'n/a's that show up in your outputs. 
#

from Bio import SeqIO 
import re,os


# Define inputs and outputs
input_file = 'results/hmmer/ef_tu/final_candidates.faa'
input_labels = 'genome_tables/Xanthomonadaceae.tsv'

accessions = {}

all_EFTus_output = 'proteins/Xanthomonadaceae/Xanthomonadaceae_all_EFTus_ver2.tsv'
uniq_EFTu22_file = 'proteins/Xanthomonadaceae/Xanthomonadaceae_uniq_elf18s_ver2.tsv'

all_EFTus = {}
uniq_EFTu22s = {}


# Create initial table with ALL EFTus in the dataset. 
# Find the EFTu22 epitope in the EFTu sequence and add it to the table. 
with open(input_labels) as strain_file:
  for line in strain_file:
    if line.startswith('GCF') or line.startswith('GCA'):
      strain_data = line.split('\t')
      accessions[strain_data[0]] = strain_data[1]

EFTu_output = open(all_EFTus_output, 'w') 

for record in SeqIO.parse(input_file, 'fasta'):
  if record.id[0:15] in accessions:
    EFTu22_seq = 'n/a'
    if re.search('[ST]IG', str(record.seq)):
      search = re.search('[ST]IG', str(record.seq))
      start = search.start()-15
      end = search.end()
      EFTu22_seq = str(record.seq)[start:end]
    EFTu_output.write(str(record.id[0:15] + '\t' + accessions[record.id[0:15]] + '\t' + record.description[16:] + '\t' + record.seq + '\t' + EFTu22_seq + '\n'))

EFTu_output.close()


# Using the previously-created table, 
# make a new table with the unique EFTu22 sequences as well as the count. 
os.system('cut -f5 ' + all_EFTus_output + ' | sort | uniq -c > proteins/Xanthomonadaceae/unique_EFTu_counts.txt')

EFTu22s_with_counts = open(uniq_EFTu22_file, 'w')

with open("proteins/Xanthomonadaceae/unique_EFTu_counts.txt") as EFTu_count_file:
  for EFTu_count in EFTu_count_file:
    with open(all_EFTus_output) as EFTu_info_file:
      for line in EFTu_info_file:
        if len(EFTu_count.strip().split(' ')) > 1:
          if EFTu_count.strip().split(' ')[1] in line:
            EFTu22s_with_counts.write(line[0:15] + '\t' + line.split('\t')[1] + '\t' + EFTu_count.strip().split(' ')[1] + '\t' + EFTu_count.strip().split(' ')[0] + '\n')
            break

EFTu22s_with_counts.close()

# os.system('rm proteins/Xanthomonadaceae/unique_EFTu_counts.txt')
        
  
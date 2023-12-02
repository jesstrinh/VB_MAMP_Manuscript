#/bin/python
#
# Given a FASTA file of MAMP candidates
# and a table of GenBank accession numbers and species names
# pull out epitope (ex: csp22) sequences from your dataset 
# and create summary tables of all sequences and unique sequences
#
# NOTE: This may not capture ALL possible csp22 sequences. 
# Manually-curate any 'n/a's that show up in your outputs. 
#

from Bio import SeqIO 
import re,os

# Define inputs and outputs
input_file = 'results/hmmer/cold_shock_prot/final_candidates.faa'
input_labels = 'genome_tables/Xanthomonadaceae.tsv'

accessions = {}

all_csps_output = 'proteins/Xanthomonadaceae/Xanthomonadaceae_all_csps_ver2.tsv'
uniq_csp22_file = 'proteins/Xanthomonadaceae/Xanthomonadaceae_uniq_csp22s_ver2.tsv'

long_csps = {}
uniq_csp22s = {}


# Create initial table with ALL CSPs in the dataset. 
# Find the csp22 epitope in the CSP sequence and add it to the table. 
with open(input_labels) as strain_file:
  for line in strain_file:
    if line.startswith('GCF') or line.startswith('GCA'):
      strain_data = line.split('\t')
      accessions[strain_data[0]] = strain_data[1]

csp_output = open(all_csps_output, 'w') 

for record in SeqIO.parse(input_file, 'fasta'):
  if record.id[0:15] in accessions:
    csp22_seq = 'n/a'
    if re.search('[KSMNPR].[YVF]G', str(record.seq)):
      search = re.search('[KSMNPR].[YVF]G', str(record.seq))
      start = search.start()-11
      end = search.end()+7
      csp22_seq = str(record.seq)[start:end]
      if len(record.seq) > 100 and str(record.seq)[start:end] not in long_csps.keys():
        long_csps[str(record.seq)[start:end]] = record.seq
    csp_output.write(str(record.id[0:15] + '\t' + accessions[record.id[0:15]] + '\t' + record.description[16:] + '\t' + record.seq + '\t' + csp22_seq + '\n'))

csp_output.close()


# Using the previously-created table, 
# make a new table with the unique csp22 sequences as well as the count. 
os.system('cut -f5 ' + all_csps_output + ' | sort | uniq -c > proteins/Xanthomonadaceae/unique_csp_counts.txt')

csp22s_with_counts = open(uniq_csp22_file, 'w')

with open("proteins/Xanthomonadaceae/unique_csp_counts.txt") as csp_count_file:
  for csp_count in csp_count_file:
    with open(all_csps_output) as csp_info_file:
      for line in csp_info_file:
        if len(csp_count.strip().split(' ')) > 1:
          if csp_count.strip().split(' ')[1] in line:
            csp22s_with_counts.write(line[0:15] + '\t' + line.split('\t')[1] + '\t' + csp_count.strip().split(' ')[1] + '\t' + csp_count.strip().split(' ')[0] + '\n')
            break

csp22s_with_counts.close()

os.system('rm proteins/Xanthomonadaceae/unique_csp_counts.txt')

long_csp_output = open('proteins/Xanthomonadaceae/csp22s_from_long_csps.txt', 'w')

for key in long_csps:
  long_csp_output.write(key + '\t' + str(long_csps[key]) + '\n')
        
long_csp_output.close()  
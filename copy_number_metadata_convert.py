#/bin/python
#
# There's inconsistent naming of the protein-annotated FASTAs and the genome FASTAs 
#
# This script is going to match the protein FASTA to the genome FASTA 
# and add a column for it 
#
# Also, add all the counts for flg, csp, and EF-Tu in one file 
# 

import os

# Input file contains the protein counts I'd like to compile
# These were created after running second_round_to_FASTA.py 
EFTu_file = 'results/hmmer/ef_tu/EFTu_counts.tsv'
csp_file = 'results/hmmer/cold_shock_prot/csp_counts.tsv'
flg_file = 'results/hmmer/flagellin/flg_counts.tsv'

temp_list = {}
temp_list_2 = {}

genome_names = 'genome_tables/genome_labels.txt'


# Output (so that you don't accidentally mess up the input file)
output_file = 'results/hmmer/all_counts_matched.tsv'


# Iterate through the input file
# When you find a filename that contains the unique GCF id,
# Add a column containing the genome filename (without the extension)
f = open(output_file, 'w')
f.write('genome\tgenome_file\tgenus\tnumber_of_flgs\tnumber_of_EFTus\tnumber_of_csps\n')

with open(flg_file) as count_tsv:
  next(count_tsv)
  for item in count_tsv:
    with open(genome_names) as name_file:
      for line in name_file:
        if line.split('\t')[0] in item:
          if int(item.split('\t')[3]) >= 5:
            temp_list[line.split('\t')[0]] = line.split('\t')[1].strip() + '\t' + item.split('\t')[2] + '\t' + '5+\n'
          else:
            temp_list[line.split('\t')[0]] = line.split('\t')[1].strip() + '\t' + item.split('\t')[2] + '\t' + item.split('\t')[3] 


# Do the same with EF-Tus
with open(EFTu_file) as count_tsv:
  next(count_tsv)
  for count in count_tsv:
    for line in temp_list:
      if line.split('\t')[0] in count:
        if int(count.split('\t')[2]) >= 5:
          temp_list_2[count.split('\t')[0]] = temp_list[count.split('\t')[0]].strip() + '\t' + '5+\n'
        else:
          temp_list_2[count.split('\t')[0]] = temp_list[count.split('\t')[0]].strip() + '\t' + count.split('\t')[2]

# Now do the same with csps
with open(csp_file) as count_tsv:
  next(count_tsv)
  for count in count_tsv:
    for line in temp_list_2: 
      if line.split('\t')[0] in count:
        if int(count.split('\t')[2]) >= 5:
          f.write(line.split('\t')[0] + '\t' + temp_list_2[count.split('\t')[0]].strip() + '\t' + '5+\n')
        else:
          f.write(line.split('\t')[0] + '\t' + temp_list_2[count.split('\t')[0]].strip() + '\t' + count.split('\t')[2])
        print(line.split('\t')[0] + '\t' + temp_list_2[count.split('\t')[0]].strip() + '\t' + count.split('\t')[2])

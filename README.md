# VB_MAMP_Manuscript
These scripts were used in my manuscript, "The perception and evolution of csp22 and elf18 from vector-borne bacterial plant pathogens".  I only uploaded custom Python scripts that I used here, but see the README for a general description of the pipeline. 


The overall goal was to pull MAMP homologs from various genome assemblies (binned by vector-borne and non-vector-borne bacterial species). After that, I applied a custom script to pull out the particular immunogenic epitope from those homologs (example: csp22 from bacterial cold shock protein) and count up a) how many MAMP copies were in each genome and b) categorize them by taxonomic class. A brief summary of the pipeline is below. 

---

MAMP homologs were pulled from assembled genomes using HMMER ver. 3.3.2 with an e-value <10^-10 and domain models corresponding to the following MAMPs:
- Bacterial flagellin - [PF00669](https://www.ebi.ac.uk/interpro/entry/pfam/PF00669/structure/PDB/) and [PF00700](https://www.ebi.ac.uk/interpro/entry/pfam/PF00700/) from the Pfam database
- Elongation factor Tu - [PF00009](https://www.ebi.ac.uk/interpro/entry/pfam/PF00009/), [PF03144](https://www.ebi.ac.uk/interpro/entry/pfam/PF03144/) and [PF03143](https://www.ebi.ac.uk/interpro/entry/pfam/PF03143/https://www.ebi.ac.uk/interpro/entry/pfam/PF03143/) from the Pfam database
- Cold shock protein - [PF00313](https://www.ebi.ac.uk/interpro/entry/pfam/PF00313/)

A list of protein IDs from the resulting table were pulled from HMMER output using shell script. This list of protein IDs was fed into the script, _first_round_to_FASTA.py_, to get a FASTA file of protein candidates for the second round of HMMER. These candidates were matched against the entire Pfam database of domain models (Pfam-A, release ver. 35.0, e-value <10^-10) to check if these hits were exclusive to MAMPs. 

Another list of protein IDs from the resulting table were pulled using shell script and fed into _second_round_to_FASTA.py_. This script not only generates a FASTA file with all your final candidates for further analysis, but will also generate summary tables on how many MAMP copies each genome has. 

We found that Liberibacters were the only vector-borne pathogen in our dataset to have a detectable copy of flagellin, so we focused on cold shock protein and elongation factor Tu instead. 

To pull out the immunogenic epitopes from each genome, we used a custom script to pull out peptide sequences based on key amino acid residues and their positions. You have to run two different scripts depending on the MAMP:
- Pulling elf18s from elongation factor Tu: _pull_ef18s.py_
- Pulling csp22s from cold shock protein: _pull_csp22s.py_
  
I manually curated the results from these scripts because not EVERY sequence of csp22 or elf18 were able to be retrieved automatically. These two scripts also generate two tables: one of csp22 epitopes grouped by taxonomic class, and one of unique epitopes with counts to later sort out epitopes by abundance in dataset. 

Lastly, since I wanted to pull all the MAMP copy number data together, I used a script to pool genome files, protein files, and copy number counts together: _copy_number_metadata_convert.py_. This was mostly to plot the data on R. 

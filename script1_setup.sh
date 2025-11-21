#!/bin/bash

#Download proteome of Acinetobacter baumannii 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/035/845/GCF_009035845.1_ASM903584v1/GCF_009035845.1_ASM903584v1_protein.faa.gz
gunzip GCF_009035845.1_ASM903584v1_protein.faa.gz

#Create new file with the  hypothetical proteins alone
awk '/^>/{p=tolower($0) ~ /hypothetical/} p' GCF_009035845.1_ASM903584v1_protein.faa > hypothetical_proteins.faa

#Filter the created file to only retain sequences longer than 50 amino acids. 
seqtk seq -L 50 hypothetical_proteins.faa > hypothetical_proteins_filtered.faa



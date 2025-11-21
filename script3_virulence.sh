#!/bin/bash

#Load required  modules 
module load blast+/2.14.1

#Set relative paths
query=$1
db=$2
out=$3

#Create output directory if necessary
mkdir -p virulence_results

#Carry out BLASTP search for virulence-related genes
blastp -query $query -db $db -out $out -evalue 1e-5 -outfmt 6 -num_threads 16

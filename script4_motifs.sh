#!/bin/bash

#Load hmmscan module
module load hmmer/3.3

#Set relative paths
input=$1
output=$2

#Carry out a protein domain search using HMMER
hmmscan --cpu 4 --domtblout $output Pfam.hmm $input

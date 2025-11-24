##======================================================================
##This script is used to process the "representative_proteins.faa" 
##imported from R Studio. It produces the "hmmscan_results.tbl" file,
##which is again imported into R Studio.
##======================================================================

#!/bin/bash

#Load required module
module load hmmer/3.3

#Set relative paths
input=$1
output=$2

#Carry out a protein domain search using HMMER
hmmscan --cpu 4 --domtblout $output Pfam.hmm $input

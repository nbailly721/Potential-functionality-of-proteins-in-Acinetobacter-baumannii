#!/bin/bash

#Load required modules
module load StdEnv/2020 
module load cd-hit/4.8.1

#Set relative paths
input=$1
output_prefix=$2
identity=$3
numbercpu=$4

#Carry out protein sequence clustering analysis using CD-HIT
cd-hit -i $input -o ${output_prefix}.faa -c $identity -T $numbercpu -M 0

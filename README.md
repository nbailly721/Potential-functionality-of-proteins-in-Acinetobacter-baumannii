                                                      ##Functional Analysis of Hypothetical Proteins in Acinetobacter baumannii##

**Description**

This project performs an end-to-end functional analysis of hypothetical proteins in Acinetobacter baumannii. It integrates protein clustering, virulence factor annotation, motif discovery using HMMER, and visualization of statistically significant motif–protein associations. The workflow processes raw protein sequences, identifies high-identity virulence hits, extracts representative proteins, analyzes motif matches, and summarizes patterns relevant to functional annotation and microbial genomics.

**Workflow Overview**
1. Data Ingestion & Preparation (R)

Load clustered proteins (clustered_proteins.faa.clstr), hypothetical proteins (hypothetical_proteins.faa), and virulence BLAST hits (virulence_results.tsv).

Inspect raw sequence metadata (headers, IDs, lengths).

Convert CD-HIT cluster output into a tidy table containing:

cluster IDs

protein IDs

sequence lengths

representative status

Clean and standardize virulence tables for merging.

2. Virulence Factor Analysis

Merge representative proteins with BLASTP virulence hits.

Filter for high-confidence results (Perc_Identity > 90%).

Summarize virulence associations for each protein:

number of hits

unique virulence factors

identity statistics

Export summary table (virulence_summary.csv).

3. Motif Discovery Preparation

Extract representative hypothetical protein sequences.

Export representative_proteins.faa for HMMER hmmscan.

4. Motif Scan Analysis

Import hmmscan_results.tbl.

Clean headers and standardize fields.

Filter for significant hits (E-value < 1e-5).

Organize motif–protein associations for visualization.

5. Visualization

Plot significant motif matches across representative proteins using ggplot2.

Use −log10(E-value) to standardize significance.

Export final figure as motif_plot.png.

**Datasets Used**

Primary Dataset: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/035/845/GCF_009035845.1_ASM903584v1/GCF_009035845.1_ASM903584v1_protein.faa.gz

Processed/Generated Files

hypothetical_proteins.faa – Filtered hypothetical protein sequences.

clustered_proteins.faa.clstr – CD-HIT protein cluster output.

`virulence_results.tsv` – BLASTP virulence factor hits.

hmmscan_results.tbl – HMMER motif scan results from representative proteins.

**Packages Used**

R Packages

tidyverse – Data manipulation & visualization

Biostrings – Sequence handling

janitor – Header cleaning

readr – File import

stringr – Pattern matching

Bash / Unix Tools

CD-HIT – Protein clustering

BLAST+ – Virulence factor annotation

HMMER – Motif/domain detection

seqtk – Sequence filtering

Note: Bash scripts assume an HPC or Linux system with module support (cd-hit/4.8.1, blast+/2.14.1, hmmer/3.3). Local users must install tools manually.

**Key Results**

motif_plot.png – Visualization of significant motif–protein relationships across representative hypothetical proteins.

virulence_summary.csv – Summary of high-identity BLASTP virulence factor hits.

**Files in This Repository**

main_analysis.R – Full R workflow (ingestion → analysis → visualization)

cluster_result.sh – CD-HIT clustering script

hypothetical_result.sh – Script to extract/clean hypothetical protein sequences

motifs_result.sh – HMMER motif scanning script

virulence_result.sh – BLASTP virulence annotation script

motif_plot.png – Motif visualization output

virulence_summary.csv – Processed virulence summary table

**Important Notes**

The workflow is fully reproducible using provided scripts and input files.

Easily adaptable to other organisms or proteomes by modifying input file paths.

Sequence and motif analysis steps are modular and can be reused independently.

**Real-World Relevance**

Supports functional annotation of hypothetical proteins in bacterial pathogens.

Identifies candidate proteins with virulence potential, aiding antimicrobial research.

Provides a reproducible computational pipeline useful for comparative proteomics.

Helps prioritize proteins for experimental validation in wet-lab studies.

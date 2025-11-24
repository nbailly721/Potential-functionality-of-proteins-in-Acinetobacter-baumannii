                                                      ##Functional Analysis of Hypothetical Proteins in Acinetobacter baumannii##

This project performs an end-to-end analysis of hypothetical proteins in Acinetobacter baumannii, integrating protein clustering, virulence factor annotation, and motif discovery. The workflow imports raw protein sequences and clustering results, identifies high-identity virulence hits, extracts representative proteins for motif scanning, processes HMMER results, and visualizes significant motif–protein associations.

🖥️ Workflow Overview

Data Ingestion & Preparation (R Script)

Load protein clusters (clustered_proteins.faa.clstr), virulence hits (virulence_results.tsv), and hypothetical protein sequences (hypothetical_proteins.faa).

Inspect raw datasets to verify headers, IDs, and sequence lengths.

Restructure protein cluster outputs into a tidy table with cluster IDs, protein IDs, sequence lengths, and representative status.

Clean and rename virulence hits for consistent merging.

Virulence Factor Analysis

Merge representative proteins with virulence hits.

Filter for high-identity hits (Perc_Identity > 90).

Summarize each protein’s virulence associations, including number of hits, unique factors, and identity statistics.

Export summary table (virulence_summary.csv) for reporting.

Motif Discovery Preparation

Extract representative protein sequences from the hypothetical protein dataset.

Export representative_proteins.faa for HMMER hmmscan analysis.

Motif Scan Analysis

Import HMMER hmmscan_results.tbl.

Clean and standardize headers.

Filter for significant motif matches (E_value < 1e-5).

Arrange and visualize protein–motif associations.

Visualization

Plot significant motif matches across representative proteins using ggplot2.

E-values are represented as -log10(E_value) for standardized significance.

Output visualization as motif_plot.png.

📁 Datasets

hypothetical_proteins.faa – Filtered hypothetical protein sequences.

clustered_proteins.faa.clstr – Protein cluster output from CD-HIT.

virulence_results.tsv – BLASTP results for virulence factor hits.

hmmscan_results.tbl – HMMER motif scan results (produced from representative proteins).

🔧 Tools & Packages

R Packages:

tidyverse – Data manipulation and visualization.

Biostrings – Protein sequence handling.

janitor – Header cleaning and standardization.

readr – Import TSV files.

stringr – String extraction and pattern matching.

Bash / Unix Tools:

CD-HIT – Protein clustering.

BLAST+ – Virulence factor annotation.

HMMER – Motif/domain scanning.

seqtk – Filtering protein sequences by length.

📊 Key Results – Figures & Tables

motif_plot.png – Significant motif matches across representative hypothetical proteins.

virulence_summary.csv – Summary table of proteins with high-identity virulence factor hits.

📂 Files in Repository

main_analysis.R – Complete R script for data ingestion, analysis, and visualization.

cluster_result.sh – Bash script to produce clustered protein output from CD-HIT.

hypothetical_result.sh – Bash script to extract and filter hypothetical proteins.

motifs_result.sh – Bash script for HMMER motif scanning.

virulence_result.sh – Bash script to run BLASTP virulence search.

motif_plot.png – Example visualization of protein–motif associations.

virulence_summary.csv – Processed summary table of virulence hits.

🧠 Notes

This workflow demonstrates a full functional analysis pipeline for bacterial hypothetical proteins.

The workflow is reproducible: raw sequences and clustering scripts are included, and outputs can be regenerated using the bash scripts and main_analysis.R.

It can be adapted to other organisms or protein datasets by updating file paths and sequence inputs.

📌 Relevance

Provides a systematic method to identify candidate proteins with virulence potential.

Highlights representative proteins for motif discovery, which can inform functional annotation studies.

Visualizes significant motif–protein relationships to support downstream experimental validation.

Can be applied to other bacterial species or genomic datasets for comparative proteomics.

🏛️ Alignment

Supports microbial genomics research, functional annotation, and computational biology reproducibility best practices.

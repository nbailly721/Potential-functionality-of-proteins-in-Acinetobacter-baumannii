## ========================================================
## Project: Functional Analysis of Hypothetical Proteins in 
## *Acinetobacter baumannii*
## Description:
##   This script analyzes hypothetical proteins in *Acinetobacter baumannii*
##   by integrating protein clustering results, virulence factor annotations,
##   and motif scan outputs. 
##   It begins by importing protein clusters, VFDB virulence hits, and 
##   hypothetical protein sequences, then restructures and merges these 
##   datasets to identify high-identity virulence associations. 
##   Representative proteins are selected for motif discovery, and HMMER 
##   hmmscan results are imported and cleaned for further analysis. 
##   Finally, the script filters significant motif matches and visualizes 
##   key motif–protein relationships.
## ========================================================

#_Set up environment -------------

install.packages('tidyverse')
install.packages('Biostrings')
install.packages('janitor')
install.packages('readr')
library('tidyverse')
library('Biostrings')
library('janitor')
library('readr')

#_Load raw data -------------

virulence_hits <- read_tsv('../data/virulence_results.tsv') 
cluster_protein <-readLines('./data/clustered_proteins.faa.clstr') 
hypothetical_proteins <- readAAStringSet('./data/hypothetical_proteins.faa') 
#Data sets imported from UNIX environment. 

#_File inspection ------------

#__Hypothetical proteins ------------

names(hypothetical_proteins) 
#To ensure that the proteins have valid IDs

width(hypothetical_proteins) 
#To determine that sequence lengths were valid

#__Virulence hits ------------

head(virulence_hits,10)
#To determine that the table had the proper headings and their respective values

#__Protein clusters ------------

head(cluster_protein,10)
#To identify what parts of the data set need to be properly formatted to simplify downstream analysis

#_Data set modification ------------

#__Protein cluster ------------

cluster_tibble <- tibble(line=cluster_protein) %>% 
  mutate(Cluster_ID = if_else(
    str_detect(line, "^>Cluster"),
    str_remove(line, "^>Cluster "),
    NA_character_)) %>%
  fill(Cluster_ID)
#Convert the data file into a tibble and create a new column "Cluster_ID" that includes the lines that start with
#"^>Cluster" and then remove that prefix. This code also propagate the Cluster_ID to all the sequences present in the cluster

cluster_final <- cluster_tibble %>%
  filter(!str_detect(line, "^>Cluster")) %>%
  mutate(
    Protein_ID = str_extract(line, "(?<=\\>)[^\\.\\s]+\\.[0-9]+"),
    Length = as.numeric(str_extract(line, "\\d+(?=aa)")),
    Representative = str_detect(line, "\\*$")
  ) %>%
  select(Cluster_ID, Protein_ID, Length, Representative)
#Modify the new data set to remove cluster headers and extract protein info: ID, length in amino acids, 
#and representative status (*). Likewise, only keep "Cluster_ID", "Protein_ID", "Length", and "Representative" columns

rm (cluster_tibble)
#Remove unnecessary pre-processed file

head(cluster_final, 10)
#Verification of correct format

#__Virulence hits ------------

virulence_hits <- virulence_hits %>%
  rename_with(~ c(
    "Protein_ID", "VF_Hit", "Perc_Identity", "Align_Length",
    "Mismatches", "Gap_Openings", "Query_Start", "Query_End",
    "Subj_Start", "Subj_End", "Evalue", "BitScore"
  ))
#Rename the headers of each column with their correct name.

head(virulence_hits,10)
#Verification of correct format

#_Analysis of virulence hits------------------

merged_virulence_cluster <- left_join(cluster_final, virulence_hits, by='Protein_ID')
#Merge the data set of the representative proteins with that of virulence hits.

merged_virulence_cluster <- merged_virulence_cluster [!is.na(merged_virulence_cluster$VF_Hit),]
#Remove the rows with NAs from the merged dataset

output_virulence_factors <- merged_virulence_cluster %>% filter (Perc_Identity > 90)
#Filter the merged data set to only keep proteins  which virulence hits had a percent identity > 90

virulence_summary <- output_virulence_factors %>%
  group_by(Protein_ID) %>%
  summarise(
    Num_VF_Hits = n (),
    VF_hits = paste(unique(VF_Hit),collapse = " ; "),
    Avg_Perc_Identity = mean(Perc_Identity),
    Max_Perc_Identity = max(Perc_Identity)
  ) 
#Summary of the the "merged_virulence_cluster" table.

print(virulence_summary)
#Verification that summary table is correctly formatted.

#_Motif functionality --------------

#__Data preparation -------------

names(hypothetical_proteins) <- str_extract(names(hypothetical_proteins), "^[^\\s]+")
#Extract/keep only the first part of the name of each protein. This is done to facilitate downstream analysis.

head(names(hypothetical_proteins))
#To ensure that the extraction step worked sucesfully

representative_hits <- cluster_final %>% 
  filter (Representative == TRUE) %>% 
  pull(Protein_ID)
#Create a new data set containing only the representative proteins. Done so that the proteins used for motif analysis are also
#representative.

match_cluster_protein <- hypothetical_proteins[names(hypothetical_proteins) %in% representative_hits]
#Create a new data set that only contains the representative proteins. Done to reduce the number of redundant proteins analyzed.


#__Analysis ----------------

class(match_cluster_protein)
names(match_cluster_protein)
match_cluster_protein[1:5]
#Ensure that format of file is correct before exporting it for analysis.

writeXStringSet(match_cluster_protein, 
                filepath = "representative_proteins.faa")
#Export file into UNIX environment (Canada Compute-Narval Cluster) for HMMER analysis. Refer to "motifs_result" script for 
#the respective Shell script

motif_results <- read_fwf(
  "../data/hmmscan_results.tbl",
  fwf_empty("../data/hmmscan_results.tbl"), skip = 3
)
print(motif_results)
#Import data set from Unix environment. Headers need to be changed for further analysis.

colnames(motif_results) <- c("target_name", "accession", "tlen", "query_name", "accession2", "qlen", "E_value", "score", 
                             "bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc", "description"
)
print(motif_results)
#Headers modified to the correct values using the original file within the UNIX environment as reference.


motif_results <- motif_results %>% 
  select(target_name,accession,query_name,E_value) %>% 
  distinct () %>% 
  mutate (E_value = as.numeric(E_value)) %>%
  filter(E_value < 1e-5) %>%
  arrange(E_value)
print(motif_results)
#Filter and organize the imported data set to retain the observations with lower e-values. Threshold chosen arbitrarily. 

#__Visualization -------------

ggplot(motif_results, aes(x = query_name, y = target_name)) +
  geom_point(aes(size = -log10(E_value)), color = "red") +
  scale_size_continuous(name = "Significance (-log10 E)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0.5, face = 'bold')
  ) +
  labs(
    x = 'Name of representative sequence',
    y = 'Type of Motif',
    title = 'Significant Motif Matches Across Hypothetical Proteins in Acinetobacter baumannii'
  )
#Graphical representation of the type of motif associated with each sequence. The e-values were plotted as a function of log(10) 
#for standardization purposes.



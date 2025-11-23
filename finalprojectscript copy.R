##***************************
## BINF6210-Assignment4-Turdus species that geographically overlap with one another have lower genetic divergence that those who do not overlap.
##
## Nicolas Bailly
##
## 2025-12-05
##
##***************************


#Do Turdus species with overlapping geographic ranges show greater genetic divergence than those in separate regions?

#_Set up environment------------------

install.packages('tidyverse')
install.packages('dplyr')
install.packages('ggplot2')
install.packages('vegan')
install.packages('Biostrings')
install.packages('leaflet')
install.packages('phytools')
install.packages('ape')
install.packages('maps')
install.packages('DECIPHER')
install.packages("phangorn")
install.packages("patchwork")
library('tidyverse')
library('dplyr')
library('Biostrings')
library('ggplot2')
library('maps')
library('leaflet')
library('phytools')
library('vegan')
library('ape')
library('DECIPHER')
library('phangorn')
library('patchwork')


#_Data insertion and visualization ------------------

turdus_raw <- read_tsv('/Users/nicolasbailly/Documents/final_project/turdus.tsv')

names(turdus_raw)
#To determine if the columns of interest of interest are present:Species name, coordinates, marker code and nucleotide sequences.

head(turdus_raw$nuc)
#To check if the nucleotide sequences are present and if they have the correct structure. This is done to prevent any downstream errors in analysis.

table(turdus_raw$marker_code)
#To determine which marker codes are present. This is done to ensure that the marker code of interest (COI-5P) is present so that it can be filtered.


#_Modification and Organization of Data ------------

turdus_coi5p <- turdus_raw %>% filter(marker_code == 'COI-5P')
#To filter out any marker code that is not COI-5P. This is done to ensure that every future analysis only uses the sub data set with that specific marker.

table(turdus_coi5p$marker_code)
#To count the number of times that the specific marker code appears in the sub data set. The high count (954) indicates that further filtering needs to be done.

turdus_coi5p_filtered <- turdus_coi5p %>%
  group_by(species) %>%
  slice_head(n=1) %>%   
  ungroup()
#Based on the code in line 42 and the manual visualization of the sub data set, some species have more than one sequence associated with the COI-5P marker. As such, this step was done to only keep one sequence per species. Such modification is necessary because downstream sequence alignment requires each species to only have one sequence.

table(turdus_coi5p_filtered$marker_code)
#This is done to ensure, through counts, that the filtering step in line 48 correctly worked. 

turdus_coi5p %>%
  filter(species == "Turdus merula") %>%
  select(species, sampleid, nuc, marker_code)
#This is done to visualize how many sequences were present for each species before the filtering steps.

turdus_coi5p_filtered %>%
  filter(species == "Turdus merula") %>%
  select(species, sampleid, nuc,marker_code)
#This is done to ensure that every different species has only one sequence after the filtering step. Done as a final check

turdus_coi5p_filtered <- turdus_coi5p_filtered %>% 
  filter(!is.na(species))
#To remove any rows that had NAs for the species value (as was row 53). The nucleotide sequence and  marker_code columns had, through visual analysis, no NAs.This step is important because the inexistence of a species name may lead to errors in sequence alignment or phylogenetic tree

unique(turdus_coi5p_filtered$species)
table(turdus_coi5p_filtered$marker_code)


#_Sequence Alignment using DECIPHER--------------------

turdus_seqs_clean <- DNAStringSet(gsub("-", "", turdus_coi5p_filtered$nuc))
#This was done to remove the "-" in the nucleotide sequences, which was giving downstream errors in alignment calculation. Hence, removed to prevent further errors.

names(turdus_seqs_clean) <- turdus_coi5p_filtered$species
#Assign a species name to the respective modified sequence. This is done to ensure that in the creation of the phylogenetic tree, the sequences match their correct species.

print(turdus_seqs_clean)
#Check to ensure that the species were correctly assigned to the sequences. This is necessary to then do the phylogenetic tree.

?AlignSeqs

alignment <- AlignSeqs(turdus_seqs_clean, iterations = 100)
#Alignment via DECIPHER using 100 iterations to increase the accuracy of the alignment. This is done to ensure that the phylogenetic tree produced dowstream is also accurate

BrowseSeqs(alignment)
#Done to visualize the aligned sequences along with the consensus nucelotide. This was also done to check that all the sequences were properly aligned.


#_Construction of phylogeny tree using Phangorn ----------

?phangorn
?phyDat
?dist.ml
?NJ
?pml
?modelTest
?optim.pml
?order

alignment_phydat <- phyDat(as.matrix(alignment), type = "DNA")
#Data set with aligned sequences converted into phyDat object. Done because the Phagorn packages requires the latter format.

distance_matrix <- dist.ml(alignment_phydat)
#Done because the production of the Neighbor-Joining tree requires pairwise distances between aligned sequences.

neighbour_tree <- NJ(distance_matrix)
#Done because the later maximum likelihood tree construction requires an initial topology or initial tree.

fit <- pml(neighbour_tree, data = alignment_phydat)
#Done because maximum likelihood tree estimation requires an initial model that combines the initial topology with the aligned sequence data.

model_test <- modelTest(fit)
##This is done to evaluate different nucleotide substitution models using AIC. The model with the lowest AIC value is considered the best fit.

print(model_test)
#To visualize more easily the AIC values

model_test_df <- as.data.frame.matrix(model_test)
#Converted to a data frame so the models can be ordered by their AIC values in the next line of code.

top_5_models <- model_test_df[order(model_test_df$AIC), ][1:5, ]
print(top_5_models) 
#To display the top 5 models ranked by lowest AIC values, helping identify the best-fitting substitution model
                
fit_gtr_model <- optim.pml(fit, model = "TIM2", optGamma = TRUE, optINV= TRUE,rearrangement = "ratchet")
#This is done to optmize the maximum likelihood tree using the best-fitting substitution and additional parameters.

plot(fit_gtr_model$tree, main = "Phylogeny of Turdus species based on COI-5p")
#This is done to visualize the inferred phylogeny within the genus and which species are more closely related.




#_Calculation of genetic divergence -------------

#__Calculation using Maximum-Likelihood Genetic Distances -----------

?dist.ml
?as.matrix

maxlikelihood_distance <- dist.ml(alignment_phydat)
#This is done to quantify pairwise genetic distances. The results are modified for visualization in the next line of code

maxlikelihood_distance <- as.data.frame(as.matrix(maxlikelihood_distance))
print(maxlikelihood_distance)
#Convert the distance object into a full matrix to visualize and use the pairwise values


print(turdus_raw$coord)

#_Geographic distribution of species ----------------

#__Modification of coordinate sub data set ---------------

turdus_coords <- turdus_raw %>% 
  select(species,coord) %>% 
  filter(!is.na(species)) %>% 
  filter(!is.na(coord)) %>%
  separate(coord, into=c("lat", "lon"), sep= ",", convert=TRUE) %>%
  mutate(
    lat = lat,
    lon = lon
  )
#This is done to remove all irrelevant columns to remove noise. Likewise, it facilitates downstream modification and analysis by removing any rows with NA values and separating the coordinates into latitude and longitude columns

turdus_coords <- turdus_coords %>%
  mutate(
    lat = as.numeric(str_replace_all(lat, "\\[|\\]|\"|\\s", "")),
    lon = as.numeric(str_replace_all(lon, "\\[|\\]|\"|\\s", ""))
  )
#This is done to ensure that the values inside the latitude and longitude cells are properly formatted, which will prevent downstream format errors.

#__Classification into overlap or no overlap category---------------

distribution_summary <- turdus_coords %>%
  group_by(species) %>%
  summarise(
    min_lat = min(lat),
    max_lat = max(lat),
    min_lon = min(lon),
    max_lon = max(lon)
  ) %>%
  ungroup()
#As there may be more than 1 observation of the same species (each with their own coordinates), this code gives an overall geographic range for every species. Such modification facilitates doing downstream pairwise comparisons where the whole geographic range of two species are compared to determine if there is an overlap.

species_pair_overlap <- expand.grid(
  species1 = as.character(distribution_summary$species),
  species2 = as.character(distribution_summary$species),
  stringsAsFactors = FALSE
) %>%
  filter(species1 != species2)
#Create all possible pairs of species (excluding self-pairs) to later compare their geographic ranges for overlap.

?paste0

species_pair_overlap <- species_pair_overlap %>%
  left_join(distribution_summary, by = c("species1" = "species")) %>%
  rename_with(~ paste0(.x, "1"), .cols = c("min_lat","max_lat","min_lon","max_lon")) %>%
  left_join(distribution_summary, by = c("species2" = "species")) %>%
  rename_with(~ paste0(.x, "2"), .cols = c("min_lat","max_lat","min_lon","max_lon"))
#This code assigns the geographic range to each species in the pair so that downstream calculations can determine which ones overlap.

species_pair_overlap <- species_pair_overlap %>%
  mutate(
    overlap = (max_lat1 >= min_lat2 & min_lat1 <= max_lat2) &  
      (max_lon1 >= min_lon2 & min_lon1 <= max_lon2)     
  )
#This is done to determine which pair of species overlap. To do so, a column called 'Overlap' was created and the repsective criteria had to be met for it to show up as 'True'.


#_Statistical test ---------------

#__Pre-analysis data set modification -------------

maxlikelihood_distance_long <- maxlikelihood_distance %>%
  rownames_to_column('species1') %>%
  pivot_longer (
      cols = -species1,
      names_to = 'species2',
      values_to = 'genetic_distance'
  ) %>%
  filter(species1 != species2)
#The data set needed to be converted into a long format in order to be correctly merged with "species_pair_overlap" data set.

maxlikelihood_distance_long_clean <- maxlikelihood_distance_long %>%
  filter(species1 < species2)
#To filter out duplicates and only remain with 1 unique pair per species combination. This is necessary to avoid inaccurate downstream statistical tests.

merged_data <- species_pair_overlap %>%
  left_join (maxlikelihood_distance_long_clean,
             by= c('species1', 'species2')) %>% 
  filter(!is.na(genetic_distance))
#This was necessary because the downstream statistical tests rely on the data being in the same data set. Removing the NAs was also necessary to prevent any errors.


#__Checking dataset consistency and missing species ------------

nrow(species_pair_overlap)
# This determines the number of rows (2070). Such number represents  all possible species pairs based on geographic ranges.
nrow(maxlikelihood_distance_long_clean)
#This determines the number of rows (1326) for which genetic distance data is available after removing duplicates and keeping only pairs with COI-5P sequences.
#The comparison between the number of rows between the two data sets to be merged suggests that some species in "species_pair_overlap" do not have available sequences. 

setdiff(unique(species_pair_overlap$species1),
        unique(maxlikelihood_distance_long_clean$species1))
setdiff(unique(species_pair_overlap$species2),
        unique(maxlikelihood_distance_long_clean$species2))
#This is done to determine that the species in the "species_pair_overlap" data set not found in the "maxlikelihood_distance_long_clean" data set is "Turdus xanthorhynchus". Hence, it serves as an extra check that the disparity in number of rows is partially due to missing species. 
# Note: The difference (744 rows) is due to: Missing sequences for some species, removal of duplicate pairs (species1 < species2), and filtering out pairs without genetic distance.

#__Determination of normality ----------------

?shapiro.test
?wilcox.test

genetic_overlap <- merged_data$genetic_distance[merged_data$overlap == TRUE]
shapiro.test(genetic_overlap) 

genetic_nooverlap <- merged_data$genetic_distance[merged_data$overlap == FALSE]
shapiro.test(genetic_nooverlap) 
#These normality tests were done to determine if the statistical tests that were going to be used to determine significant had to be non-parametric or parametric. As both tests showed a non-normal distribution, the Wilcoxon test will be used.

#__Determination of significance in difference ----------------

wilcox.test(genetic_distance ~ overlap, data = merged_data)
#This was done to determine if the difference in genetic distance between the "overlap" and "nooverlap" groups were significant. The results show that such distance is significantly different. Such step was done in order to  answer the question.

merged_data_summary <- merged_data %>%
  group_by(overlap) %>%
  summarise(
    median_distance = median(genetic_distance),
    mean_distance = mean(genetic_distance)
  )
print(merged_data_summary)
#This was done to answer the research question: The species without geographic overlap have a higher average genetic divergence..


#_Visualizations ---------------

#__Geographic overlap----------------

species_with_coi <- unique(turdus_coi5p$species)
print(species_with_coi)
#List all species with COI-5P sequences to confirm which species should be retained for mapping.

turdus_coords_filtered <- turdus_coords %>%
  filter(species %in% species_with_coi)
#Keep only coordinate records for species that have COI-5P data to ensure geographic plots match the genetic data set.

expected_species <- unique(turdus_coi5p_filtered$species)
mapped_species <- unique(turdus_coords_filtered$species)
setdiff(expected_species, mapped_species)
#This was done to determine the identify of any species  present in the "turdus_coi5p_filtered" data set, but not in the "turdus_coords_filtered" data set.

turdus_coords_missing <- turdus_coords %>%
  filter(species %in% c(
    "Turdus abyssinicus",
    "Turdus boulboul",
    "Turdus kessleri",
    "Turdus libonyana",
    "Turdus mupinensis",
    "Turdus reevei"
  ))
#These were used as an extra check to confirm that the coordinates for this species were missing. Hence, this prevents not plotting relevant species.

species_all <- unique(turdus_coords_filtered$species)
#To get a complete list of species to split into 2 graphs downstream. Such splitting is necessary because placing more than 40 species in one graph causes some of them to have the same color within the map.

n <- length(species_all)
print(n)
#Count the total number of species to determine how to separate them

group1 <- species_all[1:ceiling(n/2)]
# To plot the first graph with rouhgly half of the species

group2 <- species_all[(ceiling(n/2)+1):n]
# To create the second group with the remaining species so all species are included in one of the two plots.

coords_group1 <- turdus_coords_filtered %>% filter(species %in% group1)
coords_group2 <- turdus_coords_filtered %>% filter(species %in% group2)
#Creating these new data sets is necessary to then plot each of their respective graphs and show the geographical overlap.

?map_data
world_distribution <-map_data('world')

map1 <-ggplot() +
  geom_polygon(data = world_distribution,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey50") +
  geom_point(data = coords_group1,
             aes(x = lon, y = lat, color = species),
             size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Geographic Distribution of Turdus Species (23 Species)",
       x = "Longitude", y = "Latitude",
       color = "Species") +
  theme(plot.title=element_text(hjust=0.5, face='bold'))
#This was done to visually confirm how Turdus species are geographically distributed. Plotting the species ranges helps to determine where overlap occurs and where it does not.

map2 <- ggplot() +
  geom_polygon(data = world_distribution,
               aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey50") +
  geom_point(data = coords_group2,
             aes(x = lon, y = lat, color = species),
             size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Geographic Distribution of Turdus Species (23 Species)",
       x = "Longitude", y = "Latitude",
       color = "Species") +
  theme(plot.title=element_text(hjust=0.5, face='bold'))
#This was done to visually confirm how Turdus species are geographically distributed. Plotting the species ranges helps to determine where overlap occurs and where it does not.

combined_maps <- map1 / map2
plot(combined_maps)
#This was done to facilitate the visualization of all the geographic distributions at one.


#__Association between phylogeny and geographic overlap --------------

#___Boxplot ---------------

boxplot_overlap <- ggplot(merged_data, aes(x=overlap, y=genetic_distance, fill=overlap)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = c("FALSE" = "orange", "TRUE" = "purple")) +
  labs(
    title = 'Association between geographic overlap and genetic divergence',
    x = 'Geographic Overlap (True or False)',
    y = 'Genetic Distance',
    caption = 'The Genetic Distance is based on the Maximun Likelihood approach'
  ) +
  theme_minimal () +
  theme(
    plot.title = element_text(hjust=0.5, face = 'bold')
  )
print(boxplot_overlap)
#This was done to show how the overlap status is related to the genetic divergence within species in the group. Likewise, as previously seen, it shows how the median genetic divergence in the no-overlap group is higher than that of the overlap group.

#___Non-metric Multidimensional Scaling plot --------------

genetic_matrix <- merged_data %>%
  select(species1, species2, genetic_distance) %>%
  pivot_wider(names_from = species2, values_from = genetic_distance) %>%
  column_to_rownames("species1") %>%
  as.matrix()
#This is done to build the genetic matrix that is later used for ordination

genetic_dist <- vegdist(genetic_matrix, method = "euclidean", na.rm = TRUE)
#The distance matrix is converted into a distance object because it is needed in such format for the NMDS analysis

mds_res <- metaMDS(genetic_dist, k = 2, trymax = 100, trace = FALSE)
#The non-metric multidimensional scaling is used to convert the genetic distance into a 2D format to facilitate visualization

mds_coords <- as.data.frame(scores(mds_res))
mds_coords$species <- rownames(mds_coords)
#This is done to create a table where each species with valid genetic distance has a position in the NMDS format. Without this, the plot would not be able to plot or color the species.

mds_plot_data <- merged_data %>%
  select(species1, species2, overlap) %>%
  pivot_longer(cols = c(species1, species2), names_to = "which_species", values_to = "species") %>%
  distinct() %>%
  filter(species %in% mds_coords$species) %>%
  left_join(mds_coords, by = "species")
#This is done to integrate the data regarding the geographic overlap. In such way, the plot becomes more meaningful because it shows how genetic divergence within the genus is associated with the geographic distribution.

ggplot(mds_plot_data, aes(x = NMDS1, y = NMDS2, color = overlap)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "purple", "FALSE" = "orange")) +
  labs(title = "Association between genetic divergence and geographic overlap in Turdus species",
       x = "NMDS1", y = "NMDS2",
       color = "Geographic Overlap") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#This is done to explore the association between where species live and how genetically similar they are

nrow(mds_coords)
#Checks how many species are included in the NMDS plot, ensuring all relevant species are represented compared to other graphs. All but one were plotted.

setdiff(turdus_coords_filtered$species, mds_coords$species)
#This was done to determine the identity of the species not plotted. It was observed that the missing species is the one missing the nucleotide sequence.

#___Phylogeny-Linked Geographic Map -------------

tree <- fit_gtr_model$tree
#This was done so that the map can later show the genetic similarity between species as a measurement of genetic distances.

coords_df <- turdus_coords_filtered %>%
  group_by(species) %>%
  summarise(
    lat = mean(lat), 
    lon = mean(lon)) %>%
  ungroup() %>%
  filter(species %in% tree$tip.label) %>%
  arrange(match(species, tree$tip.label))
#This was done to give each species a single representative (average) geographical point for mapping, ensuring that the map is not overcrowded with data points. Likewise, it ensures that only species present in the phylogenetic tree are displayed in the correct order. This is necessary to accurately link geographic distribution with phylogenetic relationships.

coords_df <- coords_df %>% mutate(tree_order = seq_len(n()))
#This was done to make the color gradient represent evolutionary relationships in the map.

world_distribution <- map_data("world")
ggplot() +
  geom_polygon(data = world_distribution, aes(x = long, y = lat, group = group),
  fill = "grey90", color = "grey50") +
  geom_point(data = coords_df, aes(x = lon, y = lat, color = tree_order), 
  size = 3) +
  scale_color_viridis_c(option = "plasma", name = "Tree Tip Order") +
  labs(
  title = "Turdus Species: Phylogeny Linked to Geographic Distribution",
  x = "Longitude",
  y = "Latitude",
  caption = "Species with similar colors = More closely related genetically;
             Each point = Average lat & long location of each species;
             Lower Tree Tip Order = Species appearing lower in the tree"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#This was done to determine  if closely genetically related species are geographically clustered or sparsely distributed.


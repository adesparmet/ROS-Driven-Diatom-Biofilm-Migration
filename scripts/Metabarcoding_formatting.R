#Library packages ####
#install.packages("BiocManager")
#BiocManager::install("phyloseq", dependencies = TRUE)
library(phyloseq)
#install.packages("devtools")
#devtools::install_github("vmikk/metagMisc")
library(metagMisc)
#install.packages("iNEXT")
library(iNEXT)
#BiocManager::install("microbiome")
library(microbiome)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("ape")
library(ape)
#install.packages("here")
library(here)

# Data import ####

phyloseq_16s <- readRDS(here("data", "Metabarcoding_data", "phyloseq_16s.rds"))

phyloseq_18s <- readRDS(here("data", "Metabarcoding_data", "phyloseq_18s.rds"))

# Taxonomy table 
taxonomy_biofilm18S <- tax_table(phyloseq_18s)
taxonomy_biofilm16S <- tax_table(phyloseq_16s)

taxonomy_biofilm18S[is.na(taxonomy_biofilm18S) | taxonomy_biofilm18S == ""] <- "IncertaeSedis" #replace na by IncertaeSedis
taxonomy_biofilm16S[is.na(taxonomy_biofilm16S) | taxonomy_biofilm16S == ""] <- "IncertaeSedis"

# Unique Id
last_col18 <- ncol(taxonomy_biofilm18S)
taxonomy_biofilm18S[, last_col18] <- paste0(taxonomy_biofilm18S[, last_col18], "_n", seq_len(nrow(taxonomy_biofilm18S)))
tax_table(phyloseq_18s) <- taxonomy_biofilm18S

last_col16 <- ncol(taxonomy_biofilm16S)
taxonomy_biofilm16S[, last_col16] <- paste0(taxonomy_biofilm16S[, last_col16], "_n", seq_len(nrow(taxonomy_biofilm16S)))
tax_table(phyloseq_16s) <- taxonomy_biofilm16S

# Sample names otu_table
name_sample <-c("HL_WSB_3", "HL_WSB_2","HL_WSB_1","HL_WSB_4","HL_WSB_5","VLL_WSB_5","HP_WSB_5","VLL_WSB_2","VLL_WSB_4","VLL_WSB_1","VLL_WSB_3","HP_WSB_3","HL_SB_3","HL_SB_5","HL_SB_2","VLL_SB_2","VLL_SB_1","HL_SB_4","HL_SB_1","VLL_SB_4","VLL_SB_3","VLL_SB_5","PL_WSB_5","PL_WSB_4","PL_WSB_2","PL_WSB_3","HP_WSB_4","HP_WSB_1","PL_WSB_1","HP_WSB_2", "PL_SB_5","PL_SB_4","PL_SB_2","PL_SB_3","PL_SB_1", "HP_SB_3","HP_SB_2","HP_SB_4","HP_SB_5", "HP_SB_1")

sample_names(phyloseq_18s) <- name_sample
sample_names(phyloseq_16s) <- name_sample

# Deleting assignment errors
tax_table18 <- tax_table(phyloseq_18s)
tax_table18_filtered <- tax_table18[!(tax_table18[, "Domain"] == "Bacteria"), ]
phyloseq_18s <- prune_taxa(rownames(tax_table18_filtered), phyloseq_18s)

tax_table16 <- tax_table(phyloseq_16s)
tax_table16_filtered <- tax_table16[!(tax_table16[, "Kingdom"] == "Eukaryota"), ]
phyloseq_16s <- prune_taxa(rownames(tax_table16_filtered), phyloseq_16s)

# Rarefaction normalization ####
set.seed(12345)

biofilm18S <- phyloseq_coverage_raref(phyloseq_18s, multithread=T)
biofilm18S %>% sample_sums()
biofilm18S <- prune_taxa(taxa_sums(biofilm18S)!=0, biofilm18S) # Suppression of null lines


biofilm16S <- phyloseq_coverage_raref(phyloseq_16s, multithread=T)
biofilm16S %>% sample_sums()
biofilm16S <- prune_taxa(taxa_sums(biofilm16S)!=0, biofilm16S)


# Suppression of sequencing depth artifacts between SB and WSB ####

# Explanations : WSB is a sub-sample of SB, all WSB-specific species are rare species that appear due to the sequencing depth, so they must be removed to avoid inducing artefactual variability.

# Retrieving 18s abundances data
div_biofilm18S<- as.data.frame(otu_table(biofilm18S)) #Take biofilm18S before compositional computation
div_biofilm18S<-t(div_biofilm18S)
div_biofilm18S<-as.data.frame(div_biofilm18S)
div_biofilm18S <- cbind(Groups = sub("^[^_]+_([^_]+)_.*$", "\\1", rownames(div_biofilm18S)), div_biofilm18S) # Retrieving Groups metadata
species18_columns <- names(div_biofilm18S)[2:ncol(div_biofilm18S)]

# Sum of abundances per species for both mobilities
div_biofilm18S <- div_biofilm18S %>%
  group_by(Groups) %>%
  summarise(across(all_of(species18_columns), \(x) sum(x, na.rm = TRUE)))

# Species identification specific to WSB samples
div_biofilm18S <- div_biofilm18S %>%
  select(Groups, all_of(species18_columns)) %>%
  pivot_longer(cols = -Groups, names_to = "Species", values_to = "Abundance") %>%
  filter(Abundance != 0)

div_biofilm18S <- div_biofilm18S %>%
  group_by(Species) %>%
  summarise(
    Total_Abundance = sum(Abundance, na.rm = TRUE),
    Groups = paste(unique(Groups), collapse = "")) %>%
  arrange(desc(Total_Abundance))

table(div_biofilm18S$Groups) #result

list18_WSB <- div_biofilm18S  %>%
  filter(Groups == "WSB") %>%
  arrange(desc(Total_Abundance))
list18_WSB <- list18_WSB$Species #Retriving the list

# Deletion of species specific to WSB samples in the phyloseq object
biofilm18S@otu_table <- otu_table(biofilm18S)[!rownames(otu_table(biofilm18S)) %in% list18_WSB, ]
biofilm18S@tax_table <- tax_table(biofilm18S)[!rownames(tax_table(biofilm18S)) %in% list18_WSB, ]

final_biofilm18S <- microbiome::transform(biofilm18S, "compositional") # Transformation of abundances data between 0 and 1

# phylo_tree update 
phylo_tree18s <- rtree(n = ntaxa(biofilm18S), rooted = FALSE, tip.label = taxa_names(biofilm18S)) # Without root
phy_tree(biofilm18S) <- phylo_tree18s

      # Same filtering for 16s community

# Retrieving 16s abundances data
div_biofilm16S<- as.data.frame(otu_table(biofilm16S)) #Take biofilm16S before compositional computation
div_biofilm16S<-t(div_biofilm16S)
div_biofilm16S<-as.data.frame(div_biofilm16S)
div_biofilm16S <- cbind(Groups = sub("^[^_]+_([^_]+)_.*$", "\\1", rownames(div_biofilm16S)), div_biofilm16S) # Retrieving Groups metadata
species16_columns <- names(div_biofilm16S)[2:ncol(div_biofilm16S)]

# Sum of abundances per species for both mobilities
div_biofilm16S <- div_biofilm16S %>%
  group_by(Groups) %>%
  summarise(across(all_of(species16_columns), \(x) sum(x, na.rm = TRUE)))

# Species identification specific to WSB samples
div_biofilm16S <- div_biofilm16S %>%
  select(Groups, all_of(species16_columns)) %>%
  pivot_longer(cols = -Groups, names_to = "Species", values_to = "Abundance") %>%
  filter(Abundance != 0)

div_biofilm16S <- div_biofilm16S %>%
  group_by(Species) %>%
  summarise(
    Total_Abundance = sum(Abundance, na.rm = TRUE),
    Groups = paste(unique(Groups), collapse = "")) %>%
  arrange(desc(Total_Abundance))

table(div_biofilm16S$Groups) #result

list16_WSB <- div_biofilm16S  %>%
  filter(Groups == "WSB") %>%
  arrange(desc(Total_Abundance))
list16_WSB <- list16_WSB$Species #Retriving the list

# Deletion of species specific to WSB samples in the phyloseq object
biofilm16S@otu_table <- otu_table(biofilm16S)[!rownames(otu_table(biofilm16S)) %in% list16_WSB, ]
biofilm16S@tax_table <- tax_table(biofilm16S)[!rownames(tax_table(biofilm16S)) %in% list16_WSB, ]

final_biofilm16S <- microbiome::transform(biofilm16S, "compositional") # Transformation of abundances data between 0 and 1

# phylo_tree update 
phylo_tree16s <- rtree(n = ntaxa(biofilm16S), rooted = FALSE, tip.label = taxa_names(biofilm16S)) # Without root
phy_tree(biofilm16S) <- phylo_tree16s



# Community matrix (df_community) ####

#18s matrix
abundance18 <- as.data.frame(otu_table(final_biofilm18S))
taxa18 <- as.data.frame(tax_table(final_biofilm18S))
df_combined18 <- merge(abundance18, taxa18, by = "row.names", all.x = TRUE)

df_combined18 <- df_combined18 %>%
  mutate(Taxa = apply(select(., 42:ncol(.)), 1, function(x) paste(x, collapse = ";"))) %>%
  select(-c(42:50)) # merge taxa name
df_combined18 <- df_combined18[-1] # useless column

df_combined18 <- df_combined18 %>% select(Taxa, everything())
df_combined18 <- as.data.frame(t(df_combined18))
colnames(df_combined18) <- df_combined18[1, ] # Rownames
df_combined18 <- df_combined18[-1, ] #useless column

#16s matrix
abundance16 <- as.data.frame(otu_table(final_biofilm16S))
taxa16 <- as.data.frame(tax_table(final_biofilm16S))
df_combined16 <- merge(abundance16, taxa16, by = "row.names", all.x = TRUE)

df_combined16 <- df_combined16 %>%
  mutate(Taxa = apply(select(., 42:ncol(.)), 1, function(x) paste(x, collapse = ";"))) %>%
  select(-c(42:48)) 
df_combined16 <- df_combined16[-1] 

df_combined16 <- df_combined16 %>% select(Taxa, everything())
df_combined16 <- as.data.frame(t(df_combined16))
colnames(df_combined16) <- df_combined16[1, ] 
df_combined16 <- df_combined16[-1, ] 

#merge 18s and 16s matrix
df_community <- merge(df_combined18, df_combined16, by = "row.names", all = TRUE)

#Add metadata
df_community <- cbind(Conditions = sub("_.*", "", df_community$Row.names), df_community) #Stress conditions metadata (VLL : Very Low Light Adapted ; HL : High Light ; HP : Hydrogen Peroxide ; PL : Plamsa)
df_community <- cbind(Groups = sub("^[^_]+_([^_]+)_.*$", "\\1", df_community$Row.names), df_community) #Mobility metadata (SB : Sediment Biofilm ; WSB : Without Sediment Biofilm)
names(df_community)[names(df_community) == "Row.names"] <- "File"

#Df
df_community<- as.data.frame(df_community)
df_community[, 4:10657] <- lapply(df_community[, 4:10657], as.numeric)


# Packages
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# Load data
sponge_raw <- fread("https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00000563/pipelines/4.1/file/ERP012972_taxonomy_abundances_SSU_v4.1.tsv")

#
sponge_clean <- sponge_raw %>% 
  rename(taxonomy = `#SampleID`) %>% 
  mutate(taxonomy = str_remove_all(taxonomy, c("sk__|k__|p__|c__|o__|f__|g__|s__"))) %>%
  separate_wider_delim(cols = taxonomy, delim = ";", 
                       names = c("Domain", "Kingdom", "Phylum", 
                                 "Class", "Order", "Family", 
                                 "Genus", "Species"), 
                       too_many = "merge", too_few = "align_start") %>% 
  select(-Kingdom) %>% 
  mutate(Phylum = ifelse(Phylum == "", NA, Phylum),
         Class = ifelse(Class == "", NA, Class),
         Order = ifelse(Order == "", NA, Order),
         Family = ifelse(Family == "", NA, Family),
         Genus = ifelse(Genus == "", NA, Genus),
         Species = ifelse(Species == "", NA, Species)) %>%
  pivot_longer(cols = contains("ERR"), names_to = "Sample", values_to = "Abundance") %>% 
  filter(!Domain %in% c("Eukaryota", "Chloroplast", "Mitochondria"))

#
sponge_clean

library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0( annotations$Reference, annotations$Position, annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

lake_constraint <- fread("../data/Lake_supplementary_dataset_7.tsv") %>%
  mutate(mutation = paste0(Reference, Position, Alternate))
mdf <- merge(annotations, lake_constraint, by = "mutation")

ggplot(mdf %>% filter(syn_annotation != "other"), aes(x = syn_annotation, y = MLC_score)) +
  geom_boxplot()

mdf %>% filter(syn_annotation != "other") %>%
  filter(syn_annotation == "WCF_to_WCF")

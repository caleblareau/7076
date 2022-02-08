library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
source("01_functions.R")

# Import Phewas associations
impact <- fread("../data/cancer/data_mutations_impact.txt.gz") %>%
  filter(Variant_Type == "SNP") %>%
  mutate(mutation = ShortVariantID)

tcga <- fread("../data/cancer/data_mutations_tcga.txt.gz") %>%
  filter(Variant_Type == "SNP" & class == "somatic" & helix_freq < 0.1 & genbank_freq < 0.1) %>%
  mutate(mutation = ShortVariantID)

pcawg <- fread("../data/cancer/data_mutations_pcawg.txt.gz") %>%
  filter(Variant_Type == "SNP") %>%
  mutate(mutation = paste0(Reference_Allele, Start_Position, Tumor_Seq_Allele2))

data.frame(
  mutation = c(impact$mutation, tcga$mutation, pcawg$mutation),
  position = c(impact$Start_Position, tcga$Start_Position, pcawg$Start_Position)
) %>%
  group_by(mutation,position) %>% summarize(count_tumors = n()) -> cancer_count_df

trna_positions <- fread("../../reference_genome/output/tRNA_anticodon_positions.tsv")
trna_positions %>%
  filter(V3 %in% cancer_count_df$mutation)

merge(cancer_count_df, trna_positions, by.x = "position", by.y = "V3") %>%
  arrange((position))

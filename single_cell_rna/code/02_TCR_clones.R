library(data.table)
library(Seurat)
library(BuenColors)
library(dplyr)

# Import data
counts <- Read10X_h5("../data/T-PBMCs_Gex_29March2021_rna.rna.h5")
mtDNA_counts <- fread("../data/T-PBMCs_Gex_29March2021_7076counts.tsv")
azimuth_qc <- fread("../data/T-PBMCs_Gex_29March2021_rna.azimuth_pred.tsv")
tcrs <- fread("../data/T-PBMCs_TCR_29March2021-all_contig_annotations.csv") %>%
  filter(high_confidence == TRUE & exact_subclonotype_id == 1)

# Azimuth has already been QC'd so adjust accordingly 
mtdna_df <- mtDNA_counts %>% filter(Barcode %in% azimuth_qc$cell) %>%
  mutate(coverage = a7076 + g7076, heteroplasmy = round(g7076/(coverage + 0.00001)*100, 2))
mdf <- data.frame(merge(mtdna_df, azimuth_qc, by.x = "Barcode", by.y = "cell"))
mdf2 <- merge(mdf, tcrs, by.x = "Barcode", by.y = "barcode")

df <- mdf2 %>%
  filter(coverage > 5) %>%
  filter(heteroplasmy >99 | heteroplasmy < 1) %>%
  mutate(allele = ifelse(heteroplasmy > 99, "G7076", "A7076")) %>% 
  group_by(allele, raw_clonotype_id) %>%
  summarize(count = n())
clone_number_df <- reshape2::dcast(df, raw_clonotype_id ~ allele, fill = 0) %>% 
  arrange((G7076)) %>% data.frame()

clone_number_df_filt <- clone_number_df %>%
  mutate(prop = A7076/(A7076 + G7076)) %>%
  filter(!(prop > 0.2 & prop < 0.8))

mdf2 %>% filter(raw_clonotype_id == "clonotype4")

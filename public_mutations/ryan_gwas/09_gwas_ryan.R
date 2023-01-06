library(data.table)
library(dplyr)
library(BuenColors)
source("../code/01_functions.R")

####
# ´import nature genetics assoc
# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("M-", annotations$Position, "-",annotations$Reference, "-", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

annotations[,c("mutation","Position", "Reference", "Variant", "Symbol", "Codons", "Reference.tRNA","Variant.tRNA","syn_annotation")] %>%
  filter(!(syn_annotation == "other")) %>%
  write.table("../output/syn-table-for-ryan.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Import Phewas associations
assoc <- readxl::read_xlsx("../data/MitoPhewas_associations.xlsx", 1) %>% data.frame()
assoc$mutation <- paste0("M-", assoc$Position, "-",toupper(assoc$Allele1),"-",toupper(assoc$Allele2))
assoc_2 <- assoc
assoc_2$is_syn <- assoc_2$mutation %in% (annotations %>% filter(syn_annotation != "other") %>%
                                           pull(mutation))
assoc_2 %>%
  filter((FDR < 0.1) & is_syn) %>% pull(mutation) -> vars

assoc_2 %>% filter(FDR < 0.05 &is_syn) %>%
  filter(Meaning == "Lymphocyte count (LYMPH#)") %>%
  arrange(FDR) 

lymphvars <- c("M-11467-A-G", "M-3480-A-G","M-9698-T-C","M-11299-T-C","M-10550-A-G ")

####
df <- fread('mt_syn_binary.tsv') %>%
  filter(model == "Ryan2") %>% filter(homoplasmic_MAF < 0.05) 

df %>% filter(genotype %in% c("M-11899-T-C"))

df %>% filter(syn_annotation == "Wobble_to_Wobble") %>% head(20)
df %>% filter(syn_annotation == "Wobble_to_WCF") %>% head(20)
df %>% filter(syn_annotation == "WCF_to_Wobble") %>% head(20)

df %>% filter(pValue < 10^-6) %>% pull(Symbol) %>% table()

df %>% filter(pValue < 10^-6) %>%
  group_by(Symbol, syn_annotation, genotype) %>% 
  summarize(count = n()) %>% ungroup() %>%
  group_by( syn_annotation) %>%
  summarize(count = n()) %>% data.frame()

df %>% 
  filter(pValue < 9e-5) %>%
  mutate(Consequence = "synonymous_variant") %>% 
  process_stats_syn_type()

df %>% filter(pValue < 1e-5) %>% 
  mutate(Consequence = "synonymous_variant") %>%
  count_me()


########
dfq <- fread('mt_syn_quant.tsv') %>% data.frame()
dfq %>% filter(syn_annotation == "WCF_to_Wobble") %>% head(20)
dfq %>% filter(syn_annotation == "Wobble_to_WCF") %>% head(20)
dfq %>% filter(syn_annotation == "Wobble_to_Wobble") %>% head(20)

dfq %>% filter(geno %in% c("M-13117-A-G", "M-8655-C-T", "M-10142-C-T", "M-10188-A-G", "M-14484-T-C"))

dfq %>% filter(p < 1*10^-6) %>% 
  mutate(Consequence = "synonymous_variant") %>% 
  process_stats_syn_type()


# Try kruskal wallis test
bdf <- data.frame(
  syn_annotation = c(df$syn_annotation, dfq$syn_annotation),
  p_value = c(df$pValue, dfq$p)
)

wilcox.test(bdf$p_value[bdf$syn_annotation=="WCF_to_WCF"],
            bdf$p_value[bdf$syn_annotation!="WCF_to_WCF"])


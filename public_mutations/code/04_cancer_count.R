library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0( annotations$Reference, annotations$Position, annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

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
  mutation = c(impact$mutation, tcga$mutation, pcawg$mutation) 
) %>%
  group_by(mutation) %>% summarize(count_tumors = n()) -> cancer_count_df

# Merge the two and look at FDR < 0.05 hits (how significance has been previously defined)
cancer_merge_df <- merge(cancer_count_df, annotations, by.x = "mutation",  by.y = "mutation") %>% arrange(desc(count_tumors))
table(cancer_merge_df$syn_annotation)


# Import the GWAS data to do some combined analyses
assoc <- readxl::read_xlsx("../data/MitoPhewas_associations.xlsx", 1) %>% data.frame()
assoc$mutation <- paste0(toupper(assoc$Allele1), assoc$Position,  toupper(assoc$Allele2))

# Merge the two and look at FDR < 0.05 hits (how significance has been previously defined)
gwas_raw <- merge(assoc, annotations, by = "mutation") %>% arrange(FDR) %>% filter(FDR < 0.05)

# count only each gwas variant once
gwas <- gwas_raw %>% group_by(Consequence, Symbol, mutation, syn_annotation) %>% summarize(count = n())

# PLOT 1
# Relative frequencies of each class
cm_g <- count_me(gwas)
cm_c <- count_me(cancer_merge_df)
cm_a <- count_me(annotations)

cm_c
cm_a
cm_g

# Test statistics cancer
prop.test(cm_c$count[1], sum(cm_c$count), p = cm_a$perc[1] / 100)
prop.test(cm_c$count[2], sum(cm_c$count), p = cm_a$perc[2] / 100)
prop.test(cm_c$count[3], sum(cm_c$count), p = cm_a$perc[3] / 100)
prop.test(cm_c$count[4], sum(cm_c$count), p = cm_a$perc[4] / 100)

# Test statistics GWAS
prop.test(cm_g$count[1], sum(cm_g$count), p = cm_a$perc[1] / 100)
prop.test(cm_g$count[2], sum(cm_g$count), p = cm_a$perc[2] / 100)
prop.test(cm_g$count[3], sum(cm_g$count), p = cm_a$perc[3] / 100)
prop.test(cm_g$count[4], sum(cm_g$count), p = cm_a$perc[4] / 100)


plot1_df <- data.frame(
  variant_type = rep(cm_g$syn_annotation,3),
  value = c(cm_a$perc, cm_c$perc, cm_g$perc), 
  what = c(rep("all", 4), rep("cancer", 4), rep("gwas", 4))
)

p1 <- ggplot(plot1_df, aes(x = variant_type, y = value, fill = what)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("grey20", "firebrick", "dodgerblue4") )+
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0))+ 
  labs(x = "", y = "% of synonymous mtDNA variants") + theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../output/freqs_syn_type.pdf", width = 2.5, height = 1.8)

# Focus on WCF to Wobble 
cancer_genes <- cancer_merge_df %>% filter(syn_annotation == "WCF_to_Wobble") %>%
  group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)

all_genes <- annotations %>% filter(syn_annotation == "WCF_to_Wobble") %>%
  group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)

# group by complex
cvec <- c("CV", "CV", rep("CIV", 3), "CIII", rep("CI", 7)); names(cvec) <- all_genes$Symbol
annotations$Complex <- cvec[as.character(annotations$Symbol)]
cancer_merge_df$Complex <- cvec[as.character(cancer_merge_df$Symbol)]
gwas$Complex <- cvec[as.character(gwas$Symbol)]

# Now summarize and plot
all_complex_test <- annotations %>% filter(syn_annotation == "Wobble_to_Wobble") %>%
  group_by(Complex) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)
all_complex_test
cancer_complex_test <- cancer_merge_df %>% filter(syn_annotation == "Wobble_to_Wobble") %>%
  group_by(Complex) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)
cancer_complex_test

c(
  prop.test(cancer_complex_test$count[1], sum(cancer_complex_test$count), p = all_complex_test$perc[1] / 100)$p.value,
  prop.test(cancer_complex_test$count[2], sum(cancer_complex_test$count), p = all_complex_test$perc[2] / 100)$p.value,
  prop.test(cancer_complex_test$count[3], sum(cancer_complex_test$count), p = all_complex_test$perc[3] / 100)$p.value,
  prop.test(cancer_complex_test$count[4], sum(cancer_complex_test$count), p = all_complex_test$perc[4] / 100)$p.value
)


# Gwas complex association
gwas_wcf_to_wobble <- gwas %>% filter(syn_annotation == "WCF_to_Wobble") %>%
  group_by(Complex) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)

prop.test(gwas_wcf_to_wobble$count[1], sum(gwas_wcf_to_wobble$count), p = all_complex_wcf_to_wobble$perc[1] / 100)
prop.test(gwas_wcf_to_wobble$count[2], sum(gwas_wcf_to_wobble$count), p = all_complex_wcf_to_wobble$perc[2] / 100)
prop.test(gwas_wcf_to_wobble$count[3], sum(gwas_wcf_to_wobble$count), p = all_complex_wcf_to_wobble$perc[3] / 100)
prop.test(gwas_wcf_to_wobble$count[4], sum(gwas_wcf_to_wobble$count), p = all_complex_wcf_to_wobble$perc[4] / 100)



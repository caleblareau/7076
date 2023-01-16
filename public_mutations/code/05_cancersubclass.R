library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0( annotations$Reference, annotations$Position, annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)
anno_aa_df <- annotations %>% filter(syn_annotation == "WCF_to_Wobble") %>% 
  group_by(Amino.acids) %>% summarize(count = n()) %>% mutate(pct = count / sum(count)*100) 


# Import cancer associations
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
  cancer = c(impact$maintype, tcga$maintype, pcawg$project)
) -> cancer_df

table(cancer_df$cancer)

# Merge the two and look at FDR < 0.05 hits (how significance has been previously defined)
cancer_merge_df <- merge(cancer_df, annotations, by.x = "mutation",  by.y = "mutation")
table(cancer_merge_df$syn_annotation)

cm_a <- count_me(annotations)
# Function to count things
process_stats_syn <- function(df){
  cm_a <- count_me(annotations)
  cm_t <- count_me(df %>% mutate(syn_annotation = factor(syn_annotation, levels = cm_a$syn_annotation)) )
  cm_t$all_perc <- (cm_a$perc)
  cm_t$pvalue <- c(prop.test(cm_t$count[1], sum(cm_t$count), p = cm_a$perc[1] / 100)$p.value,
                   prop.test(cm_t$count[2], sum(cm_t$count), p = cm_a$perc[2] / 100)$p.value,
                   prop.test(cm_t$count[3], sum(cm_t$count), p = cm_a$perc[3] / 100)$p.value,
                   prop.test(cm_t$count[4], sum(cm_t$count), p = cm_a$perc[4] / 100)$p.value)
  cm_t
  
}
process_stats_syn(cancer_merge_df)

# Function to count things
count_simple_tissue <- function(df, tissue_keep){
  print(tissue_keep)
  cm_a <- count_me(annotations)
  cm_t <- count_me(df %>% mutate(syn_annotation = factor(syn_annotation, levels = cm_a$syn_annotation)) %>% filter(cancer == tissue_keep))
  data.frame(tissue_keep, WCF_to_Wobble = as.numeric(cm_t[2,2]), Wobble_to_WCF = as.numeric(cm_t[3,2]),
             total = sum(cm_t[,2]))
}

# PLOT 1
# Relative frequencies of each class
tissue_count <- sort(table(cancer_merge_df$cancer))

lapply(names(tissue_count)[as.numeric(tissue_count) > 20], function(one_tissue){
  cm_g <- count_simple_tissue(cancer_merge_df,one_tissue)
  cm_g %>%
    mutate(pvalue = prop.test(Wobble_to_WCF, total, p = (22.0 / 100))$p.value)
}) %>% rbindlist() %>% data.frame() %>%  mutate(ratio = (Wobble_to_WCF + 0.1)/(total + 0.1)) %>%
  filter(total > 20) %>%
  arrange(desc(ratio))-> all_assoc
all_assoc
all_assoc$padj <- p.adjust(all_assoc$pvalue)
all_assoc %>% arrange(padj, desc(ratio)) %>%
  mutate(rank = 1:n()) %>%
  ggplot(aes(x = rank, y = -log10(padj), color = ratio*100)) + labs(x = "Rank sorted tissues", color = "") + 
  geom_point() + scale_color_gradientn(colours = jdb_palette("solar_extra")) +
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") ->p0
cowplot::ggsave2(p0, file = "../output/gtex_tissues.pdf", width = 1.8, height = 1.8)

library(ggbeeswarm)
p1bar <- ggplot(all_df_count, aes(x = syn_annotation, y = perc, fill = syn_annotation, shape = cohort)) +
  geom_bar(data = cm_c, stat = "identity", position = "dodge", color = "black",) +
  geom_quasirandom() + scale_shape_manual(values=c(0, 1, 2,3))+
  scale_fill_manual(values = jdb_palette("corona")[1:4] )+
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), limits = c(0,45))+ 
  labs(x = "", y = "% of synonymous mtDNA variants") + theme(legend.position = "none")
p1bar
cowplot::ggsave2(p1bar, file = "../output/freqs_syn_type.pdf", width = 2, height = 1.5)

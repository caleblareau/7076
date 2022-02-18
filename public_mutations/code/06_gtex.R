library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0(annotations$Position, "_",annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

# Import gtex mutations
gtex <- readRDS("../data/GTEx_allTissueMutationsDonors_N1.rds")

# Liftover bp
hg19_to_RCS <- function(bp){
  bp[bp %in% c(310, 312, 16184)] <- 0
  bpo <- ifelse(bp < 309, bp,
                ifelse(bp == 311, 310,
                       ifelse(bp < 3108, bp -2,
                              ifelse(bp < 16183, bp -1, bp -2))))
}
sdf <- stringr::str_split_fixed(gtex$Mutation, "_", 2)
gtex$RCS_pos <- hg19_to_RCS(as.numeric(as.character(sdf[,1])))
gtex$mutation <- paste0(gtex$RCS_pos, "_", sdf[,2])

m_gtex_df <- merge(gtex, annotations, by = "mutation")

# Function to count things
process_stats_syn_type <- function(df, tissue_keep){
  print(tissue_keep)
  cm_a <- count_me(annotations)
  cm_t <- count_me(df %>% mutate(syn_annotation = factor(syn_annotation, levels = cm_a$syn_annotation)) %>% filter(tissue == tissue_keep))
  cm_t$all_perc <- (cm_a$perc)
  cm_t$pvalue <- c(prop.test(cm_t$count[1], sum(cm_t$count), p = cm_a$perc[1] / 100)$p.value,
                   prop.test(cm_t$count[2], sum(cm_t$count), p = cm_a$perc[2] / 100)$p.value,
                   prop.test(cm_t$count[3], sum(cm_t$count), p = cm_a$perc[3] / 100)$p.value,
                   prop.test(cm_t$count[4], sum(cm_t$count), p = cm_a$perc[4] / 100)$p.value)
  cm_t
  
}
# PLOT 1
# Relative frequencies of each class
gtex_tissue_count <- sort(table(m_gtex_df$tissue))

lapply(names(gtex_tissue_count)[as.numeric(gtex_tissue_count) > 20], function(one_tissue){
  cm_g <- process_stats_syn_type(m_gtex_df,one_tissue)
  cm_g$tissue <- one_tissue
  cm_g
}) %>% rbindlist() %>% data.frame() %>% arrange(pvalue) -> all_assoc
all_assoc



library(ggbeeswarm)
p1bar <- ggplot(all_df_count, aes(x = syn_annotation, y = perc, fill = syn_annotation, shape = cohort)) +
  geom_bar(data = cm_c, stat = "identity", position = "dodge", color = "black",) +
  geom_quasirandom() + scale_shape_manual(values=c(0, 1, 2,3))+
  scale_fill_manual(values = jdb_palette("corona")[1:4] )+
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), limits = c(0,45))+ 
  labs(x = "", y = "% of synonymous mtDNA variants") + theme(legend.position = "none")
p1bar
cowplot::ggsave2(p1bar, file = "../output/freqs_syn_type.pdf", width = 2, height = 1.5)

# Focus on WCF to Wobble 
cancer_genes <- cancer_merge_df %>% filter(syn_annotation == "WCF_to_Wobble") %>%
  group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)

all_genes <- annotations %>% filter(syn_annotation == "WCF_to_Wobble") %>%
  group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)

# group by complex
cvec <- c("CV", "CV", rep("CIV", 3), "CIII", rep("CI", 7)); names(cvec) <- all_genes$Symbol
annotations$Complex <- cvec[as.character(annotations$Symbol)]
gwas$Complex <- cvec[as.character(gwas$Symbol)]

process_complex_association <- function(variant_set, what_mutations = "WCF_to_Wobble"){
  
  if(what_mutations == "stop_gained"){
    annotations$syn_annotation <- annotations$Consequence
    cancer_merge_df2$syn_annotation <- cancer_merge_df2$Consequence
  }
  
  data.frame(
    mutation = c(variant_set) 
  ) %>%
    group_by(mutation) %>% summarize(count_tumors = n()) -> cancer_count_df2
  cancer_merge_df2 <- merge(cancer_count_df2, annotations, by.x = "mutation",  by.y = "mutation") %>% arrange(desc(count_tumors))
  
  
  cancer_merge_df2$Complex <- cvec[as.character(cancer_merge_df2$Symbol)]
  
  # Now summarize and plot
  all_complex_test <- annotations %>% filter(syn_annotation == what_mutations) %>%
    group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)
  all_complex_test
  cancer_complex_test <- cancer_merge_df2 %>% mutate(Symbol = as.factor(Symbol)) %>% filter(syn_annotation == what_mutations) %>%
    group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100) %>%
    tidyr::complete(Symbol, fill = list(perc = 0, count = 0))
  cancer_complex_test
  
  #pval_vec <- sapply(1:13, function(i){
  #  prop.test(cancer_complex_test$count[i], sum(cancer_complex_test$count), p = all_complex_test$perc[i] / 100)$p.value
  #})
  names(all_complex_test) <- c("Symbola", "counta", "perca")
  cbind(cancer_complex_test, all_complex_test)
}

# count number of associatios per gene
process_complex_association <- function(variant_set, what_mutations = "WCF_to_Wobble"){
  
  data.frame(
    mutation = c(variant_set) 
  ) %>%
    group_by(mutation) %>% summarize(count_tumors = n()) -> cancer_count_df2
  cancer_merge_df2 <- merge(cancer_count_df2, annotations, by.x = "mutation",  by.y = "mutation") %>% arrange(desc(count_tumors))
  
  if(what_mutations == "stop_gained"){
    annotations$syn_annotation <- annotations$Consequence
    cancer_merge_df2$syn_annotation <- cancer_merge_df2$Consequence
  }
  
  if(what_mutations == "all_others"){
    annotations$syn_annotation <- ifelse(annotations$syn_annotation %in% c("WCF_to_Wobble", "other"), "a", "all_others")
    cancer_merge_df2$syn_annotation <- ifelse(cancer_merge_df2$syn_annotation %in% c("WCF_to_Wobble", "other"), "a", "all_others")
  }
  
  
  # Now summarize and plot
  all_complex_test <- annotations %>% filter(syn_annotation == what_mutations) %>%
    group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100)
  all_complex_test
  cancer_complex_test <- cancer_merge_df2 %>% mutate(Symbol = as.factor(Symbol)) %>% filter(syn_annotation == what_mutations) %>%
    group_by(Symbol) %>% summarize(count = n()) %>% mutate(perc = count / sum(count)*100) %>%
    tidyr::complete(Symbol, fill = list(perc = 0, count = 0))
  cancer_complex_test
  
  #pval_vec <- sapply(1:13, function(i){
  #  prop.test(cancer_complex_test$count[i], sum(cancer_complex_test$count), p = all_complex_test$perc[i] / 100)$p.value
  #})
  names(all_complex_test) <- c("Symbola", "counta", "perca")
  cbind(cancer_complex_test, all_complex_test)
}

syn <- process_complex_association(c(impact$mutation,tcga$mutation, pcawg$mutation), "WCF_to_Wobble")

cbind(
  process_complex_association(c(impact$mutation), "WCF_to_Wobble"),
  process_complex_association(c(tcga$mutation), "WCF_to_Wobble"),
  process_complex_association(c(pcawg$mutation), "WCF_to_Wobble")
)

missense <- process_complex_association(c(impact$mutation,tcga$mutation, pcawg$mutation), "stop_gained")

# Plot the two 
plot_df <- data.frame(
  symbol = syn$Symbol, 
  rat_syn = (syn$perc-syn$perca), 
  rat_missense = (missense$perc-missense$perca),
  syn = syn$count/syn$counta,
  missense = missense$count/missense$counta
)

library(ggrepel)
p1 <- ggplot(plot_df , aes(x = (rat_syn), y = (rat_missense), label = symbol)) +
  geom_smooth(method=lm, se=FALSE) +
  geom_point(size = 0.5) + 
  labs(x = "WCF->Wobble (%Obs. - %Exp.)", y = "Missense Variants(%Obs. - %Exp.)") +
  pretty_plot(fontsize = 7) + L_border() 
p1
cowplot::ggsave2(p1, file = "../output/per_gene_mutation_rates.pdf", width = 1.5, height = 1.5)

cor.test((syn$perc-syn$perca), (missense$perc-missense$perca), method = "spearman")
cor.test((syn$count/syn$counta), (missense$count/missense$counta), method = "spearman")

# Now do some control tests
syn_wobble_to_wcf <- process_complex_association(c(impact$mutation,tcga$mutation, pcawg$mutation), "all_others")
missense <- process_complex_association(c(impact$mutation,tcga$mutation, pcawg$mutation), "other")

plot_df2 <- data.frame(
  symbol = syn_wobble_to_wcf$Symbol, 
  rat_syn = (syn_wobble_to_wcf$perc-syn_wobble_to_wcf$perca), 
  rat_missense = (missense$perc-missense$perca),
  syn = syn_wobble_to_wcf$count / syn_wobble_to_wcf/syn$counta,
  missense = missense$count/missense$counta
)
cor.test((syn_wobble_to_wcf$perc-syn_wobble_to_wcf$perca), (missense$perc-missense$perca), method = "spearman")
cor.test((syn_wobble_to_wcf$count/syn_wobble_to_wcf$counta), (missense$count/missense$counta), method = "spearman")

p2 <- ggplot(plot_df2, aes(x = (rat_syn), y = (rat_missense), label = symbol)) +
  geom_smooth(method=lm, se=FALSE) +
  geom_point(size = 0.5) + 
  labs(x = "Wobble->WCF (%Obs. - %Exp.)", y = "Missense Variants(%Obs. - %Exp.)") +
  pretty_plot(fontsize = 7) + L_border() 

p2
cowplot::ggsave2(p2, file = "../output/per_gene_mutation_rates-supplement.pdf", width = 1.5, height = 1.5)

gwas_gene <- process_complex_association(gwas$mutation, "Wobble_to_WCF")
cor.test((gwas_gene$perc-gwas_gene$perca), (missense$perc-missense$perca), method = "spearman")




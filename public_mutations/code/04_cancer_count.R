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


process_stats_syn_type2 <- function(df){
  cm_t <- count_me(df)
  cm_t
  
}

# Define functions to slice and dice variant class
filter_for_Gbois <- function(df){
  df %>% filter(Amino.acids %in% c("I", "F","C", "S")) %>%
    filter(!(Amino.acids == "S") | grepl("GCU", Reference.tRNA)) %>% #refine serine; either matches another AA or has the right anticodon
    count_me
}

filter_for_tau <- function(df){
  df %>% filter(Amino.acids %in% c("L", "W", "E", "K", "Q")) %>%
    filter(!(Amino.acids == "L") | grepl("UAA", Reference.tRNA)) %>% #refine leucine; either matches another AA or has the right anticodon
    count_me
}

filter_for_super_wobble <- function(df){
  df %>% filter(Amino.acids %in% c("G", "R", "A", "T", "P","V", "S", "L")) %>%
    filter(!(Amino.acids == "S") | grepl("UGA", Reference.tRNA)) %>% #refine serine; either matches another AA or has the right anticodon
    filter(!(Amino.acids == "L") | grepl("UAG", Reference.tRNA)) %>% #refine lysine; either matches another AA or has the right anticodon
    
    count_me
}


filter_for_Q34 <- function(df){
  df %>% filter(Amino.acids %in% c("Y", "H", "N", "D")) %>%
    count_me
}

classes_df <- rbind(
  filter_for_super_wobble(cancer_merge_df) %>% mutate(class = "uracil", what = "obs."),
  filter_for_super_wobble(annotations) %>% mutate(class = "uracil", what = "exp."),
  
  filter_for_Gbois(cancer_merge_df) %>% mutate(class = "guanine", what = "obs."),
  filter_for_Gbois(annotations) %>% mutate(class = "guanine", what = "exp."),
  
  filter_for_Q34(cancer_merge_df) %>% mutate(class = "queuosine", what = "obs."),
  filter_for_Q34(annotations) %>% mutate(class = "queuosine", what = "exp."),
  
  filter_for_tau(cancer_merge_df) %>% mutate(class = "taurine-uracil", what = "obs."),
  filter_for_tau(annotations) %>% mutate(class = "taurine-uracil", what = "exp.")
)

p1 <- ggplot(classes_df, aes(y = perc, x = what, fill = syn_annotation)) + 
  geom_bar(stat = "identity", color = "black") + facet_grid(~class) +
  scale_fill_manual(values = jdb_palette("corona")[1:4] ) +
  pretty_plot(fontsize = 7) +
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.position = "none") + labs(x = "", y = '% of synonymous mutations')
cowplot::ggsave2(p1, file = "../output/split_classes.pdf", width = 4.2, height = 1.8)

# do statistical tests 
get_pvalue <- function(class_x){
  classes_df %>%
    filter(class == class_x) %>%
    data.frame %>% 
    reshape2::dcast(.,formula = syn_annotation ~ what, value.var = "count") %>%
    dplyr::select(observed,expected) %>% data.matrix() %>% chisq.test() %>% str
}
get_pvalue ("guanine")
get_pvalue ("uracil")
get_pvalue ("queuosine")
get_pvalue ("taurine-uracil")



#######



# Import the GWAS data to do some combined analyses
assoc <- readxl::read_xlsx("../data/MitoPhewas_associations.xlsx", 1) %>% data.frame()
assoc$mutation <- paste0(toupper(assoc$Allele1), assoc$Position,  toupper(assoc$Allele2))

# Merge the two and look at FDR < 0.05 hits (how significance has been previously defined)
gwas_raw <- merge(assoc, annotations, by = "mutation") %>% arrange(FDR) %>% filter(FDR < 0.05)

# count only each gwas variant once
gwas <- gwas_raw %>% group_by(Consequence, Symbol, mutation, syn_annotation) %>% summarize(count = n())

make_cancer_df <- function(variant_set){
  data.frame(
    mutation = c(variant_set) 
  ) %>%
    group_by(mutation) %>% summarize(count_tumors = n()) -> cancer_count_df2
  cancer_merge_df2 <- merge(cancer_count_df2, annotations, by.x = "mutation",  by.y = "mutation") %>% arrange(desc(count_tumors))
  cancer_merge_df2
}

process_stats_syn_type <- function(df){
  cm_a <- count_me(annotations)
  cm_t <- count_me(df)
  cm_t$all_perc <- cm_a$perc
  cm_t$pvalue <- c(prop.test(cm_t$count[1], sum(cm_t$count), p = cm_a$perc[1] / 100)$p.value,
                   prop.test(cm_t$count[2], sum(cm_t$count), p = cm_a$perc[2] / 100)$p.value,
                   prop.test(cm_t$count[3], sum(cm_t$count), p = cm_a$perc[3] / 100)$p.value,
                   prop.test(cm_t$count[4], sum(cm_t$count), p = cm_a$perc[4] / 100)$p.value)
  cm_t
  
}
# PLOT 1
# Relative frequencies of each class
cm_g <- process_stats_syn_type(gwas)
cm_c <- process_stats_syn_type(make_cancer_df(c(impact$mutation, tcga$mutation, pcawg$mutation)))

cm_c$cohort <- ""

# Now compute them individually
impact_df_count <- process_stats_syn_type(make_cancer_df(c(impact$mutation))); impact_df_count$cohort <- "impact"
tcga_df_count <- process_stats_syn_type(make_cancer_df(c(tcga$mutation))); tcga_df_count$cohort <- "tcga"
pcawg_df_count <- process_stats_syn_type(make_cancer_df(c(pcawg$mutation))); pcawg_df_count$cohort <- "pcawg"
all_df_count <- rbind(impact_df_count, tcga_df_count, pcawg_df_count)
all_df_count

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




library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
library(VariantAnnotation)
source("01_functions.R")
vcf <- readVcf( "../../reference_genome/data/gnomad.genomes.v3.1.sites.chrM.vcf.bgz", verbose = FALSE )

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$id <- paste0(annotations$Reference,annotations$Position, annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

simple_gnomad <- data.frame(rowRanges(vcf)@ranges,
                            id=vcf@info$variant_collapsed,
                            n_homo = vcf@info$AC_hom, n_het = vcf@info$AC_het) %>% filter(width == 1)

# Now merge
mdf <- merge(simple_gnomad, annotations, by = "id")
gnomad_het <- mdf %>% filter(n_het >=1)
gnomad_homo <- mdf %>% filter(n_homo >=1)
gnomad_homo_rare <- mdf %>% filter(n_homo >=1 & n_homo < 100)
gnomad_homo_common <- mdf %>% filter( n_homo > 100)

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

process_stats_syn_type2 <- function(df){
  cm_t <- count_me(df)
  cm_t
  
}

process_stats_syn_type(gnomad_homo_rare)
process_stats_syn_type(gnomad_homo_common)


# Make it continuous
count_me_gnomad_cutoff <- function(df, cutoff = 1){
  print(cutoff)
  df %>% 
    filter(n_homo >= cutoff) %>% 
    filter(Consequence == "synonymous_variant") %>%
    group_by(syn_annotation) %>%
    summarize(count = n()) %>%
    tidyr::complete(syn_annotation, fill = list(count = 0)) %>%
    mutate(perc = count/sum(count)*100) %>%
    filter(syn_annotation == "Wobble_to_Wobble") %>%
    pull(perc)
}

lapply(c(1:10, 15, 20, 30, 40, 50, 75, 100, 150, 200,  300, 400, 500), function(threshold){
  data.frame(n_carriers = threshold, propWobbleToWCF = count_me_gnomad_cutoff(mdf, threshold))
}) %>% rbindlist() %>%
  ggplot(aes(x = n_carriers, y = propWobbleToWCF)) + 
  geom_hline(yintercept = 22, linetype = 2) +
  geom_line(color = "green3") + 
  scale_x_log10() + 
  geom_point() +
  scale_y_continuous(limits = c(15,45)) +
  pretty_plot() + L_border() +
  labs(x = "Number of carriers of homoplasmic variant",
       y = "% Wobble --> WCF")



# Make it continuous
count_me_gnomad_cutoff_het <- function(df, cutoff = 1){
  print(cutoff)
  df %>% 
    filter(n_het >= cutoff) %>% 
    filter(Consequence == "synonymous_variant") %>%
    group_by(syn_annotation) %>%
    summarize(count = n()) %>%
    tidyr::complete(syn_annotation, fill = list(count = 0)) %>%
    mutate(perc = count/sum(count)*100) %>%
    filter(syn_annotation == "Wobble_to_WCF") %>%
    pull(perc)
}

lapply(c(1:5), function(threshold){
  data.frame(n_carriers = threshold, propWobbleToWCF = count_me_gnomad_cutoff_het(mdf, threshold))
}) %>% rbindlist() %>%
  ggplot(aes(x = n_carriers, y = propWobbleToWCF)) + 
  geom_hline(yintercept = 22, linetype = 2) +
  geom_line(color = "green3") + 
  scale_x_log10() + 
  geom_point() +
  scale_y_continuous(limits = c(15,45)) +
  pretty_plot() + L_border() +
  labs(x = "Number of carriers of homoplasmic variant",
       y = "% Wobble --> WCF")


library(readxl)
library(BuenColors)
source("01_functions.R")

# Import Jacob's work now for gtex
annotations2 <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations2$mutation <- paste0(annotations2$Position, "_",annotations2$Variant)
annotations2$syn_annotation <- annotate_synonomous(annotations2)

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

m_gtex_df <- merge(gtex, annotations2, by = "mutation")

process_stats_syn_type(m_gtex_df)

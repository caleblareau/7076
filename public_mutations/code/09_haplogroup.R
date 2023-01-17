library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
library(VariantAnnotation)
source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$id <- paste0(annotations$Position, annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

# Data from here:
# https://github.com/seppinho/haplogrep-cmd/tree/master/src/main/resources
# Last column is number of genomes
haplo <- fread("../data/haplogrep.data.txt")  %>% 
  filter(V6 >= 10)
haplo$id <- gsub("!", "", haplo$V1)
minimal_df <- data.frame(id = unique(haplo$id))

# Now merge
mdf <- merge(minimal_df, annotations, by = "id") %>%
  unique()

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

mdf %>% filter(syn_annotation != "other") %>% dim()
process_stats_syn_type(mdf)


write.table(mdf %>% filter(syn_annotation != "other") %>% arrange(desc(syn_annotation)), 
            file = "../output/SuppTable_Haplogroup_annotation.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


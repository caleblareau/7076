library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)

source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0( annotations$Position, annotations$Reference, ">", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

# Import Phewas associations
assoc <- readxl::read_xlsx("../data/20210504_Suppl_Table_2_mtDNA_per_CLL.xlsx", 1) %>% data.frame()

make_cancer_df <- function(vec, unique = TRUE){
  mdf <- merge(data.frame(mutation = vec), annotations, by.x = "mutation",  by.y = "mutation")
  if(unique){
    mdf %>% group_by(mutation) %>% summarize(count_tumors = n()) -> cancer_count_df 
  } else{
    
  }
}

penter <- make_cancer_df(assoc[[1]])
penter_merge_df_all <- merge(penter, annotations, by.x = "mutation",  by.y = "mutation") %>% arrange(desc(count_tumors))

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

process_stats_syn_type(penter_merge_df_all)
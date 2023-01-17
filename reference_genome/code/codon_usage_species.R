library(data.table)
library(dplyr)
library(stringr)
library(BuenColors)

# Set up plotting preferences
order_df <- data.frame(
  rep(rep(c("T", "C", "A", "G"), each = 4),4),
  rep(c("T", "C", "A", "G"), each = 16),
  rep(c("T", "C", "A", "G"), 16))

collapse_df <- function(order_df){
  
  codon_order <- sapply(1:dim(order_df)[1], function(i){
    paste0(order_df[i,1], order_df[i,2], order_df[i,3])
  })
  codon_order
}
codon_order <- collapse_df(order_df)
codon_annotation_mito <- fread("../data/reference_tRNA_anticodons-wStop.tsv")
codon_annotation_nuc <- fread("../data/NUCLEAR-reference_tRNA_anticodons-wStop.txt")

wcf_codons <- codon_annotation_mito %>% filter(AA.matchbool == "*") %>%
  pull(Codon.dna)

species <- gsub("data_", "", gsub("_mitochondrial_codon_statistics.tsv", "", list.files("../data/mito_codonstatsdb/")))

lapply(species, function(x){
  dt <- fread(paste0("../data/mito_codonstatsdb/data_",x,"_mitochondrial_codon_statistics.tsv")) %>% data.frame()
  total <- sum(dt[["Total.num"]])
  wcf <- sum(dt[dt[["CODON"]] %in% wcf_codons,][["Total.num"]])
  data.frame(x, wcf, total) %>%
    mutate(prop = wcf/total*100)
}) %>%  rbindlist() %>% arrange(prop) %>%
  filter(x %in% c(9606, 10090))

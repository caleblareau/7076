library(dplyr)
library(BuenColors)
library(ggforce)
library(data.table)
library(stringr)

# Set up coordinates
df <- data.frame(pos = 1:16299, angle = seq(0, (2*pi), by = (2*pi)/(16299-1)))

# Import Jacob's work
annotations <- fread("../reference/mouse_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("m", annotations$Position, annotations$Ref, ">", annotations$Alt)
annotations$syn_annotation <- annotations$classify
table(annotations$syn_annotation)

# Look at specific amino acid changes
annotations %>% filter(syn_annotation == "WCF_to_WCF") %>% pull(Amino_acids) %>% table()
annotations %>% filter(syn_annotation == "Wobble_to_Wobble") %>% pull(Amino_acids) %>% table()
annotations %>% filter(syn_annotation == "WCF_to_Wobble") %>% pull(Amino_acids) %>% table()
annotations %>% filter(syn_annotation == "Wobble_to_WCF") %>% pull(Amino_acids) %>% table()

pos_count_df <- annotations %>% group_by(Position, syn_annotation) %>% summarize(count = n()) %>%
  filter(syn_annotation != "other")
mmdf <- merge(df, pos_count_df, by.x = "pos", by.y = 'Position')
mmdf$syn_annotation_int <- as.integer(as.factor(mmdf$syn_annotation))

library(ggbeeswarm)
p1 <- ggplot(mmdf, aes(x = syn_annotation_int, y = pos, color = syn_annotation))+
  geom_quasirandom(size = 0.1) + coord_polar(theta = "y") +
  theme_void() + 
  xlim(-5,5) + ylim(1, 16299) +
  scale_color_manual(values = jdb_palette("corona"))

#cowplot::ggsave2(p1, file = "../output/possible_mutations_mouse.pdf", width = 3, height = 3)

# creat lookup vec
lookup_vec <- annotations$syn_annotation; names(lookup_vec) <- annotations$mutation

######
all_muts <- fread("../data/12863_2021_1005_MOESM5_ESM.csv")
all_muts$heteroplasmy <- all_muts$mut_allele_count/(all_muts$mut_allele_count+all_muts$ref_allele_count)
table(all_muts$tissue)
all_muts$mutation <- paste0("m", all_muts$ref_num, all_muts$ref, ">", all_muts$mutated_base)
all_muts$syn_annotation <- lookup_vec[as.character(all_muts$mutation)]
all_muts <- all_muts[!is.na(all_muts$syn_annotation),]
dim(all_muts)

# enumerate
count_me <- function(df){
  df %>% 
    group_by(syn_annotation) %>%
    summarize(count = n()) %>%
    tidyr::complete(syn_annotation, fill = list(count = 0)) %>%
    mutate(perc = count/sum(count)*100)
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
process_stats_syn_type(annotations)
process_stats_syn_type(all_muts %>% filter(som_germ == "Somatic"))
process_stats_syn_type(all_muts %>% filter(tissue == "Liver")%>% filter(heteroplasmy > 0.01) %>% filter(som_germ == "Somatic"))
process_stats_syn_type(all_muts %>% filter(heteroplasmy > 0.01) %>% filter(som_germ == "Germline"))

lapply(unique(all_muts$old_animal), function(x){
  print(x)
  process_stats_syn_type(all_muts %>% filter(old_animal == x))
  
})

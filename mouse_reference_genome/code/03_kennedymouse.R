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


# create lookup vec
lookup_vec <- annotations$syn_annotation; names(lookup_vec) <- annotations$mutation

######
all_muts <- fread("../data/Sanchez_Contreras_etal_2022-Supplementary-File-2.csv")
all_muts <- all_muts %>%
  filter(ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T"))

all_muts$mutation <- paste0("m", all_muts$start, all_muts$ref, ">", all_muts$alt)
all_muts$syn_annotation <- lookup_vec[as.character(all_muts$mutation)]
all_muts <- all_muts[!is.na(all_muts$syn_annotation),]
dim(all_muts)

uniq_muts <- data.frame(
  mutation = unique(all_muts$mutation)
)
uniq_muts$syn_annotation <- lookup_vec[as.character(uniq_muts$mutation)]

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
process_stats_syn_type(all_muts )
process_stats_syn_type(uniq_muts )

#####
all_muts$tissue <- str_split_fixed(all_muts$sample, "_", 5)[,3]
all_muts$treatment <- str_split_fixed(all_muts$sample, "_", 5)[,2]
all_muts$library <- str_split_fixed(all_muts$sample, "_", 5)[,1]

# Import annotaitons
anno_df <- fread("../data/kennedy_annotation.csv")
m45 <- anno_df %>% filter(Age == 4.5) %>% pull(MouseID)
m26 <- anno_df %>% filter(Age == 26) %>% pull(MouseID)

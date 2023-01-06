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
process_stats_syn(m_gtex_df)

# Function to count things
count_simple_tissue <- function(df, tissue_keep){
  print(tissue_keep)
  cm_a <- count_me(annotations)
  cm_t <- count_me(df %>% mutate(syn_annotation = factor(syn_annotation, levels = cm_a$syn_annotation)) %>% filter(tissue == tissue_keep))
  data.frame(tissue_keep, WCF_to_Wobble = as.numeric(cm_t[2,2]), Wobble_to_WCF = as.numeric(cm_t[3,2]),
             total = sum(cm_t[,2]))
}

count_me(m_gtex_df %>% mutate(syn_annotation = factor(syn_annotation, levels = cm_a$syn_annotation)) %>% filter(tissue == "Muscle - Skeletal"))

# PLOT 1
# Relative frequencies of each class
gtex_tissue_count <- sort(table(m_gtex_df$tissue))

lapply(names(gtex_tissue_count)[as.numeric(gtex_tissue_count) > 20], function(one_tissue){
  cm_g <- count_simple_tissue(m_gtex_df,one_tissue)
  cm_g %>%
    mutate(pvalue = prop.test(Wobble_to_WCF, total, p = (22.0 / 100))$p.value)
}) %>% rbindlist() %>% data.frame() %>%  mutate(ratio = (Wobble_to_WCF + 0.1)/(total + 0.1)) %>%
  filter(total > 20) %>%
  arrange(desc(ratio))-> all_assoc
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

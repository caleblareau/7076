library(data.table)
library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(ggbeeswarm)
source('https://gist.githubusercontent.com/caleblareau/baee9629b9bf4c8ada1a833174ddef3e/raw/7e280c170128404789e0e62a1e1ef0dce1bdb09b/variant_calling.R')
source('../../public_mutations/code/01_functions.R')

anno <- fread("../../public_mutations/data/functional_variant_tRNA_anticodon_table.tsv") %>%
  mutate(mutation = paste0(Position, Reference, ">", Variant))
anno$syn <- anno %>% annotate_synonomous()

rb <- data.matrix(read.table("riboseqbedgraph.tsv", header = TRUE))
cpm <- round(t(t(rb)/colSums(rb)*1000000))
celltype_vec <- c(rep("DermalFibroblasts", 3), rep("Myoblasts", 4), rep("Myocytes", 2), rep("HelaS3",6), rep("HEK293", 3))
colnames(cpm) <- make.unique(rep(celltype_vec, each = 2))
cpm <- cpm[,c(FALSE, TRUE)]
colnames(rb)[c(FALSE, TRUE)]

# Determine mutations specific to each via RNA-seq
se <- readRDS("rna5.rds")

assays(call_mutations_mgatk(se))[["allele_frequency"]] %>%
  data.matrix() -> dm
assays(call_mutations_mgatk(se))[["coverage"]] %>%
  data.matrix() -> mrna_coverage
mrna_coverage_cpm <- round(t(t(mrna_coverage)/colSums(mrna_coverage)*1000000 *3)) # multiply by 3 b/c each position is represented 3 times
colnames(mrna_coverage_cpm) <- paste0("cov", colnames(mrna_coverage_cpm))
boo_keep <- rowSums(dm>0.99) >= 1 & rowSums(dm<0.01) >= 1 
mdf <- merge(data.frame(round(dm[boo_keep , ] )), 
             anno[,c("mutation", "syn", "Symbol","Position")],
             by.x = "row.names", by.y = "mutation") %>% filter(syn != "other")

# Identify candidate mutations
f1 <- mdf %>% filter(HelaS3!=HEK293)
subset_cpm_ribo <- cpm[f1 %>% pull(Position),celltype_vec %in% c("HEK293", "HelaS3")]
comp_df <- data.frame(f1, mean=round(rowMeans(subset_cpm_ribo),0), subset_cpm_ribo) %>%
  arrange(mean)

plot_diff_df <- comp_df %>% arrange(desc(mean)) %>% mutate(rank = 1:n())
ggplot(plot_diff_df, aes(rank, y = mean/100, color = syn)) +
  geom_point() +
  scale_color_manual(values = jdb_palette("corona")[1:4] ) +
  scale_x_continuous(breaks = c(1,5,10)) +
  pretty_plot(fontsize = 6) + L_border() +
  theme(legend.position = "none") +
  labs(x = "Synonymous differences between Hela and HEK293", y = "MitoRibo-seq cp10k coverage") -> pA
cowplot::ggsave2(pA, file = "all_muts.pdf", width = 1.7, height = 1.7)

comp_df <- comp_df %>% filter(mean > 200)
sapply(1:dim(comp_df), function(idx){
  wilcox.test(as.numeric(comp_df[idx,grepl("HelaS3", colnames(comp_df))])[-1],
              as.numeric(comp_df[idx,grepl("HEK293", colnames(comp_df))])[-1])$p.value # - 1 is to get rid of genotype column 
}) %>% round(4)


# Look at RNA coverage for a sense of translational efficiency; looks ok
mrna_coverage_cpm[c( "12372G>A", "9540T>C"),]

# Now build data frames for statistical tests and plotting
# transform from CPM to cp10k for consistency
data.frame(
  cellline= c(rep("HEK293", 3), rep("Hela", 6)),
  value = c( c(113,135,132), c(207,188,135,216,335,522))/100) -> df_12372

df_12372 %>%
  ggplot(aes(x = cellline, y = value)) +
  geom_boxplot(width = 0.7) +
  geom_quasirandom() + ggtitle("m.12372G>A") +
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "MitoRibo-seq cp10k coverage") -> pA
cowplot::ggsave2(pA, file = "m12372GA.pdf", width = 1.2, height = 1.7)
wilcox.test(value ~ cellline, df_12372)

data.frame(
  cellline= c(rep("HEK293", 3), rep("Hela", 6)),
  value = c( c(183,122,114),c(145,224,76,262,455,334))/100) -> df_9540
wilcox.test(value ~ cellline, df_9540)

df_9540 %>%
  ggplot(aes(x = cellline, y = value)) +
  geom_boxplot(width = 0.7) +
  geom_quasirandom() + ggtitle("m.9540T>C") + 
  pretty_plot(fontsize = 7) + L_border() +
  labs(x = "", y = "MitoRibo-seq cp10k coverage") -> pB
cowplot::ggsave2(pB, file = "m9540TC.pdf", width = 1.2, height = 1.7)



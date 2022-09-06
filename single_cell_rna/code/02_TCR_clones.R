library(data.table)
library(Seurat)
library(BuenColors)
library(dplyr)

# Import data
counts <- Read10X_h5("../data/T-PBMCs_Gex_29March2021_rna.rna.h5")
mtDNA_counts <- fread("../data/T-PBMCs_Gex_29March2021_7076counts.tsv")
azimuth_qc <- fread("../data/T-PBMCs_Gex_29March2021_rna.azimuth_pred.tsv")
tcrs <- fread("../data/T-PBMCs_TCR_29March2021-all_contig_annotations.csv") %>%
  filter(high_confidence == TRUE & exact_subclonotype_id == 1)

# Azimuth has already been QC'd so adjust accordingly 
mtdna_df <- mtDNA_counts %>% filter(Barcode %in% azimuth_qc$cell) %>%
  mutate(coverage = a7076 + g7076, heteroplasmy = round(g7076/(coverage + 0.00001)*100, 2))
mdf <- data.frame(merge(mtdna_df, azimuth_qc, by.x = "Barcode", by.y = "cell"))
mdf2 <- merge(mdf, tcrs, by.x = "Barcode", by.y = "barcode")

df <- mdf2 %>%
  filter(coverage > 5) %>%
  filter(predicted.celltype.l1 %in% c("CD4 T", "CD8 T")) %>%
  mutate(ct = substr(predicted.celltype.l1, 1, 3)) %>%
  filter(heteroplasmy >99 | heteroplasmy < 1) %>%
  mutate(allele = ifelse(heteroplasmy > 99, "G7076", "A7076")) %>% 
  group_by(allele, raw_clonotype_id,ct) %>%
  summarize(count = n())
table(df$ct)

# Clone df analysis
clone_number_df <- reshape2::dcast(df[,c("allele", "count", "raw_clonotype_id", "ct")],
                                   raw_clonotype_id ~ allele + ct, value.var = "count", fill = 0) %>% 
  arrange((A7076_CD8)) %>% data.frame()


# Look at all clones
clone_df <- reshape2::melt(clone_number_df, id.vars = c("raw_clonotype_id")) %>%
  filter(value > 0) %>%
  group_by(variable) %>%
  summarize(mean_clone_size = mean(value), sem_clone_size = sqrt(var(value)/n()), n(), median(value))

clone_df$variable <- factor(as.character(clone_df$variable), c("A7076_CD4", "G7076_CD4", "A7076_CD8", "G7076_CD8"))
p1 <- clone_df %>%
  ggplot(aes(x = variable, y = mean_clone_size)) +
  geom_bar(color = "black", stat = "identity", fill = "lightgrey") +
  geom_errorbar(aes(ymin=mean_clone_size-sem_clone_size, ymax=mean_clone_size+sem_clone_size), width=.2,
                position=position_dodge(.9)) +
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Mean clone size")
p1
cowplot::ggsave2(p1, file = "../plots/TCR_exansion.pdf", width = 2.5, height = 2)

# do statistical test
reshape2::melt(clone_number_df, id.vars = c("raw_clonotype_id")) %>%
  filter(value > 0) %>%
  filter(variable %in% c("A7076_CD8", "G7076_CD8")) %>%
  t.test(data = ., value ~ variable)

fdf8 <- reshape2::melt(clone_number_df, id.vars = c("raw_clonotype_id")) %>%
  filter(value > 1)  %>%
  filter(variable %in% c("G7076_CD8", "A7076_CD8"))

fdf4 <- reshape2::melt(clone_number_df, id.vars = c("raw_clonotype_id")) %>%
  filter(value > 1)  %>%
  filter(variable %in% c("G7076_CD4", "A7076_CD4"))


wilcox.test(value ~ variable, data = fdf4)
wilcox.test(value ~ variable, data = fdf8)


ggplot(fdf4, aes(x = value, color = variable)) +
  stat_ecdf() + scale_x_log10()

ggplot(fdf8, aes(x = value, color = variable)) +
  stat_ecdf() + scale_x_log10()

ks.test(fdf4$variable, fdf4$value)
ks.test(fdf8$variable, fdf8$value)


# Look at the top clone / celltype
datamatrixcount <- data.matrix(clone_number_df[,2:5])
mc <- max.col(datamatrixcount)
max_clone_df <- data.frame(
  raw_clonotype_id = clone_number_df$raw_clonotype_id[mc],
  celltype = colnames(datamatrixcount)[mc],
  number = matrixStats::rowMaxs(datamatrixcount)
)

max_clone_df %>%
  filter(number > 3) %>%
  group_by(celltype) %>%
  summarize(mean_clone_size = mean(number), sem_clone_size = sqrt(var(number)/n()), n()) 


#
clone_number_df_filt <- clone_number_df %>%
  mutate(prop = A7076/(A7076 + G7076)) %>%
  filter(!(prop > 0.2 & prop < 0.8)) %>%
  filter(A7076+G7076 > 3)

clone_number_df_filt$max_allele <- ifelse(clone_number_df_filt$A7076 > clone_number_df_filt$G7076, "A7076", "G7076")
clone_number_df_filt$allele_count <- ifelse(clone_number_df_filt$A7076 > clone_number_df_filt$G7076,
                                            clone_number_df_filt$A7076, clone_number_df_filt$G7076)
library(ggbeeswarm)
ggplot(clone_number_df_filt, aes(x = max_allele, y = allele_count)) +
  geom_violin() + 
  scale_y_log10()

mdf2 %>% filter(raw_clonotype_id == "clonotype4")

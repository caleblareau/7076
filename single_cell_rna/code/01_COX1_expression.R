library(data.table)
library(Seurat)
library(BuenColors)
library(dplyr)

# Import data
counts <- Read10X_h5("../data/T-PBMCs_Gex_29March2021_rna.rna.h5")
mtDNA_counts <- fread("../data/T-PBMCs_Gex_29March2021_7076counts.tsv")
azimuth_qc <- fread("../data/T-PBMCs_Gex_29March2021_rna.azimuth_pred.tsv")

# Azimuth has already been QC'd so adjust accordingly 
mtdna_df <- mtDNA_counts %>% dplyr::filter(Barcode %in% azimuth_qc$cell) %>%
  mutate(coverage = a7076 + g7076, heteroplasmy = round(g7076/(coverage + 0.00001)*100, 2))
mdf <- data.frame(merge(mtdna_df, azimuth_qc, by.x = "Barcode", by.y = "cell"))
rownames(mdf) <- mdf$Barcode
counts_filtered <- counts[,mdf$Barcode]

# Now Do Seurat things
so <- CreateSeuratObject(counts = counts_filtered, project = "ansu", meta.data = mdf)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- ScaleData(so)
so@meta.data$COX1_normalized <- so@assays$RNA@data["MT-CO1",]
so@meta.data$COX1_raw <- so@assays$RNA@counts["MT-CO1",]

# Add an oxphos module score

# Get scores
oxphos_pathway <- c("SNORD138","MIR4691","COX17","MIR7113","ATP5PD","ATP5MG","UQCR11","COX4I1",
                    "COX5B","COX6A1","COX6A2","COX6B1","COX6C","COX7A1","COX7A2","COX7B","COX7C",
                    "COX8A","COX11","COX15","UQCRQ","DMAC2L","SLC25A4","SLC25A5","SLC25A6","UQCR10",
                    "NDUFS7","NDUFA1","NDUFA2","NDUFA3",
                    "NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFAB1","NDUFB1",
                    "NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                    "NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFV1","NDUFS4","NDUFS5","NDUFS6",
                    "NDUFS8","NDUFV2","NDUFV3","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
                    "ATP5PB","ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5PF","ATP5PO","NDUFA12",
                    "SCO1","SDHA","SDHB","SDHC","SDHD","SURF1","UCP1","UCP2","UCP3","UQCRB","UQCRC1",
                    "UQCRC2","UQCRFS1","UQCRH","SLC25A14","COX7A2L","COX5A","ATP5IF1","SLC25A27","ATP5MF")
so <- AddModuleScore(so, list(oxphos_pathway),name = c("oxphos"))

# Just do things with a data frame for simplicity
df <- so@meta.data
count_vec <- sort(table(df$predicted.celltype.l2))
keep_celltypes <- names(count_vec[count_vec > 250])

# Filter and modify the data frame for downstream use
filt_df <- df %>% dplyr::filter(coverage >= 2) %>%
  dplyr::filter(predicted.celltype.l2 %in% keep_celltypes) %>%
  mutate(assign = case_when(
    heteroplasmy > 95 ~ "7076G", 
    heteroplasmy < 5 ~ "7076A",
    TRUE ~ "heteroplasmic"
  ))

p0 <- ggplot(filt_df %>% dplyr::filter(assign != "heteroplasmic") %>%
               dplyr::filter(predicted.celltype.l2 %in% c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TEM")), 
             aes(x = assign, y = log2(COX1_raw + 1))) + 
  facet_wrap(~predicted.celltype.l2, nrow = 1) + 
  geom_violin(outlier.shape = NA) + pretty_plot(fontsize = 8) +
  pretty_plot() + labs(x = "7076 variant status", y = "COX1 expression (log2UMIs + 1)")

#cowplot::ggsave2(p0, file = "../plots/COX1-expression.pdf",  width = 4, height = 1.8)

pOXPHOS <- ggplot(filt_df %>% dplyr::filter(assign != "heteroplasmic") %>%
               dplyr::filter(predicted.celltype.l2 %in% c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TEM")), 
             aes(x = assign, y = oxphos1)) + 
  facet_wrap(~predicted.celltype.l2, nrow = 1) + 
  geom_violin(outlier.shape = NA) + pretty_plot(fontsize = 8) +
  pretty_plot() + labs(x = "7076 variant status", y = "OXPHOS module score")
#cowplot::ggsave2(pOXPHOS, file = "../plots/pOXPHOS-expression.pdf",  width = 4, height = 1.8)


ctdf <- data.frame(
  celltype = unique(filt_df$predicted.celltype.l2)
)
ctdf$pvalue <- sapply(ctdf$celltype, function(ct){
  ssdf <- filt_df %>%
    dplyr::filter(predicted.celltype.l2 %in% ct) %>% dplyr::filter(assign != "heteroplasmic")
  wilcox.test(log2(ssdf$COX1_raw + 1) ~ ssdf$assign)$p.value
})
ctdf$pvalue

p1 <- filt_df %>%
  group_by(predicted.celltype.l2) %>%
  summarize(mean7076a_g = round(mean(heteroplasmy),2), cell_count = n()) %>%
  arrange(desc(mean7076a_g)) %>%
  mutate(rank = 1:n())  %>%
  ggplot(aes(x = predicted.celltype.l2, y = mean7076a_g, size = cell_count)) +
  geom_point() + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = "mean 7076A>G%") + scale_y_continuous(limits = c(24, 60)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_size_continuous(range = c(0.3,2.5))

p1

#cowplot::ggsave2(p1, file = "../plots/Heteroplasmy.pdf",  width = 3, height = 1.8)

filt_df %>%
  group_by(predicted.celltype.l2) %>%
  summarize(mean7076a_g = round(mean(heteroplasmy),2), cell_count = n()) %>%
  arrange(desc(mean7076a_g)) %>%
  write.table(sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#-------------
# Look at DEGs in the data

# Now do a loop to look at DE genes
lapply(keep_celltypes, function(ct){
  print(ct)
  # Get only those cells for a given celltype 
  filt_df %>%
    dplyr::filter(predicted.celltype.l2 == ct) %>%
    dplyr::filter(assign %in% c("7076A", "7076G")) -> fdf
  sos <- CreateSeuratObject(counts = counts[,rownames(fdf)], project = "ansu", meta.data = fdf)
  sos <- NormalizeData(sos, normalization.method = "LogNormalize", scale.factor = 10000)
  fm <- FindMarkers(sos, ident.1 = "7076G", ident.2 = "7076A", group.by = "assign", logfc.threshold = 0.1)
  fm$celltype <- ct
  fm$gene <- rownames(fm)
  fm
}) %>% rbindlist() %>% data.frame() -> all_diff_df

all_diff_df[!grepl("^TR", all_diff_df$gene),] %>%
  dplyr::filter(celltype %in% c("CD8 TEM", "CD8 Naive")) %>% arrange((p_val))

all_diff_df[!grepl("^TR", all_diff_df$gene),] %>%
  dplyr::filter(p_val_adj < 0.1)

# look at top genes that aren't TCR genes
pV <- all_diff_df[!grepl("^TR", all_diff_df$gene),]

all_diff_df[!grepl("^TR", all_diff_df$gene),] %>% 
  dplyr::filter(celltype %in% c("CD8 TEM", "CD8 Naive")) %>%
  arrange(p_val_adj) %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = celltype, label = gene)) +
  pretty_plot(fontsize = 8) + theme(legend.position = "none") +
  geom_text() + facet_wrap(~celltype) +
  scale_color_manual(values = jdb_palette("corona")[c(4,1)])


#cowplot::ggsave2(pV, file = "../plots/DEvolcanoes_exansion.pdf", width = 4.1, height = 2)



# Look at all celltypes
filt_df_cd8tem <- filt_df %>% dplyr::filter(predicted.celltype.l2 == 'CD8 TEM')
sos <- CreateSeuratObject(counts = counts[,rownames(filt_df_cd8tem)], project = "ansu", meta.data = filt_df_cd8tem)
sos <- NormalizeData(sos, normalization.method = "LogNormalize", scale.factor = 10000)
fm <- FindMarkers(sos, ident.1 = "7076A", ident.2 = "7076G", group.by = "assign", logfc.threshold = 0.1)
head(fm)
mean(sos@assays$RNA@counts["TIGIT",filt_df_cd8tem$assign == "7076A"])
mean(sos@assays$RNA@counts["TIGIT",filt_df_cd8tem$assign == "7076G"])

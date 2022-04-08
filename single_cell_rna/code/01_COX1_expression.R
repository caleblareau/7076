library(data.table)
library(Seurat)
library(BuenColors)
library(dplyr)

# Import data
counts <- Read10X_h5("../data/AnsuPBMCs_Gex_29March2021_rna.rna.h5")
mtDNA_counts <- fread("../data/AnsuPBMCs_Gex_29March2021_7076counts.tsv")
azimuth_qc <- fread("../data/AnsuPBMCs_Gex_29March2021_rna.azimuth_pred.tsv")

# Azimuth has already been QC'd so adjust accordingly 
mtdna_df <- mtDNA_counts %>% filter(Barcode %in% azimuth_qc$cell) %>%
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

# Just do things with a data frame for simplicity
df <- so@meta.data
count_vec <- sort(table(df$predicted.celltype.l2))
keep_celltypes <- names(count_vec[count_vec > 250])

# Filter and modify the data frame for downstream use
filt_df <- df %>% filter(coverage >= 2) %>%
  filter(predicted.celltype.l2 %in% keep_celltypes) %>%
  mutate(assign = case_when(
    heteroplasmy > 95 ~ "7076G", 
    heteroplasmy < 5 ~ "7076A",
    TRUE ~ "heteroplasmic"
  ))

ggplot(filt_df, aes(x = assign, y = log2(COX1_raw + 1))) + 
  facet_wrap(~predicted.celltype.l2) + 
  geom_boxplot(outlier.shape = NA) +
  pretty_plot() + labs(x = "7076 variant status", y = "COX1 expression (log2 + 1 expression)")

filt_df %>%
  group_by(predicted.celltype.l2) %>%
  summarize(mean7076a_g = round(mean(heteroplasmy),2), cell_count = n()) %>%
  write.table(sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#-------------
# Look at DEGs in the data

# Now do a loop to look at DE genes
lapply(keep_celltypes, function(ct){
  print(ct)
  # Get only those cells for a given celltype 
  filt_df %>%
    filter(predicted.celltype.l2 == ct) %>%
    filter(assign %in% c("7076A", "7076G")) -> fdf
  sos <- CreateSeuratObject(counts = counts[,rownames(fdf)], project = "ansu", meta.data = fdf)
  sos <- NormalizeData(sos, normalization.method = "LogNormalize", scale.factor = 10000)
  fm <- FindMarkers(sos, ident.1 = "7076A", ident.2 = "7076G", group.by = "assign", logfc.threshold = 0.1)
  fm$celltype <- ct
  fm$gene <- rownames(fm)
  fm
}) %>% rbindlist() %>% data.frame() -> all_diff_df

# look at top genes
all_diff_df %>% arrange(p_val_adj)

# Look at all celltypes
filt_df_cd8tem <- filt_df %>% filter(predicted.celltype.l2 == 'CD8 TEM')
sos <- CreateSeuratObject(counts = counts[,rownames(filt_df_cd8tem)], project = "ansu", meta.data = filt_df_cd8tem)
sos <- NormalizeData(sos, normalization.method = "LogNormalize", scale.factor = 10000)
fm <- FindMarkers(sos, ident.1 = "7076A", ident.2 = "7076G", group.by = "assign", logfc.threshold = 0.1)
head(fm)
mean(sos@assays$RNA@counts["TIGIT",filt_df_cd8tem$assign == "7076A"])
mean(sos@assays$RNA@counts["TIGIT",filt_df_cd8tem$assign == "7076G"])

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
so <- FindVariableFeatures(so)
so <- RunPCA(so)
so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.5)
so <- RunUMAP(so, dims = 1:20)

so@meta.data$COX1_normalized <- so@assays$RNA@data["MT-CO1",]
so@meta.data$COX1_raw <- so@assays$RNA@counts["MT-CO1",]

# Add an oxphos module score
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
so$m7076 <- so@meta.data %>%mutate(assign = case_when(
  coverage < 2 ~ "low_coverage",
  heteroplasmy > 95 ~ "m7076G", 
  heteroplasmy < 5 ~ "m7076A",
  TRUE ~ "heteroplasmic"
)) %>% pull(assign)

DimPlot(so, label = TRUE, group.by = c("seurat_clusters"), split.by = "m7076")
FeaturePlot(so, c("KLRD1", "GZMK", "SELL", "TNFRSF18", "FCGR3A"))

FindMarkers(so, "5", "2")

so@meta.data %>% group_by(m7076, seurat_clusters) %>% summarize(count = n()) %>%
  filter(m7076 %in% c("m7076A", "m7076G")) %>% ungroup() %>% reshape2::dcast(seurat_clusters~m7076, fill = "count") %>%
  data.frame() %>%
  mutate(log2_GA = round(log2(as.numeric(m7076G)/as.numeric(m7076A)), 2))

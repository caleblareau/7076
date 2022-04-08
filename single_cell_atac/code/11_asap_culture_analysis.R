library(Seurat)
library(Signac)
library(BuenColors)
library(dplyr)
library(data.table)
library(Matrix)

sc <- fread("../data/asap_culture_data/TEM_culture-ATAC_hg38_v20-mtMask.singlecell.csv")
mito <- readRDS("../data/asap_culture_data/TEM_culture-ATAC_hg38_v20-7076.rds")
h5 <- Read10X_h5("../data/asap_culture_data/TEM_culture-ATAC_hg38_v20.peakbycell.h5")

mito <- mito %>% filter(cov7076 >= 10 & cov7076 < 200)
mdf <- merge(mito, sc, by.x = "cell", by.y = "barcode")

# Import ADT
import_kite_counts_stim <- function(){
  mtx <- fread(paste0("../data/asap_culture_data/tag-data/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/asap_culture_data/tag-data/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]], "-1")
  colnames(matx) <- paste0(fread(paste0("../data/asap_culture_data/tag-data/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}

tags <- import_kite_counts_stim()
rat_control <- colSums(tags[grep("control", rownames(tags)),])/colSums(tags)
tags <- tags[,colSums(tags) >= 200 & rat_control < 0.05]
mdf <- mdf %>% filter(cell %in% colnames(tags))

# Now create the Seurat object
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- 'UCSC'

# Do the Seurat object now
chrom_assay <- CreateChromatinAssay(
  counts = h5[grepl("^chr", rownames(h5)),mdf$cell],
  sep = c(":", "-"),
  #genome = 'hg38',
  fragments = '../../../m7076-large-data-files/fragments-files/Ansu_TEM_culture-ATAC_hg38_v20.fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
)

rownames(mdf) <- mdf$cell
so <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = mdf
)
rm(chrom_assay)
rm(h5)

# Add ADT data
so[["ADT"]] <- CreateAssayObject(counts = tags[,rownames(mdf)])
clr_function <- function(x) {
  return(log1p(x = x/(exp(x = sum(log1p(x = x[x > 
                                                0]), na.rm = TRUE)/length(x = x)))))
}

so@assays$ADT@data = t(data.matrix(clr_function(t(so@assays$ADT@data))))
DefaultAssay(so) <- "ADT"
so <- so %>% FindVariableFeatures(nfeatures = 50, selection.method = "mvp", nbin = 3) 
VariableFeatures(so) %>% sort()

# Do ADT analyses
so <- so%>% ScaleData() %>% RunPCA(reduction.name = 'apca')
rm(tags)

# Now do chromatin data
DefaultAssay(so) <- "peaks"
so <- RunTFIDF(so) %>%
  FindTopFeatures( min.cutoff = 'q50') %>%
  RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:20, reduction.name = "peaks.umap")
so <- so %>% FindNeighbors( reduction = 'lsi', dims = 2:20, graph.name = 'peaks_lsi')
so <- so %>%  FindClusters( graph.name = "peaks_lsi", resolution = 0.8)

DimPlot(object = so, label = TRUE, reduction = "peaks.umap") 

library(viridis)
FeaturePlot(so, features = "af7076", reduction = "peaks.umap") +
  scale_color_viridis()

#so <- FindMultiModalNeighbors(
#  so, reduction.list = list("lsi", "apca"), 
#  dims.list = list(2:20, 1:10), modality.weight.name = "weight"
#)

# Now do weighted nearest neighbors
#so <- RunUMAP(so, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#so <- FindClusters(so, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)

DimPlot(object = so, label = TRUE, reduction = "wnn.umap") 

DefaultAssay(so) <- "ADT"
FeaturePlot(so, features = c("CD4", "CD8A", "CD45RA", "CD45RO", "ITGB1", "HLA-DR","KLRK1", "KLRG1", "af7076"), max.cutoff = "q95", 
            reduction = "peaks.umap") &
  scale_color_viridis()

data.frame(
  cor = cor(t(data.matrix(so@assays$ADT@data)), so$af7076)) %>%
  arrange(desc(cor))

# Annotate
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
Annotation(so) <- annotations

summary(so@meta.data$af7076)

# Identify markers associated with clusters
lapply(1:dim(clrd)[1], function(i){
  data.frame(stat = kruskal.test(clrd[i,], as.factor(so@meta.data$seurat_clusters))$statistic,
             marker = rownames(clrd)[i])
}) %>% rbindlist() %>% arrange(desc(stat)) -> marker_assoc

head(marker_assoc, 20)




DimPlot(so, label = TRUE)

so@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(prop_zero_cells = round(mean(af7076 < 0.01), 3)) %>% 
  write.table(row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# View gene scores
gene.activities <- GeneActivity(so)
so[['RNA']] <- CreateAssayObject(counts = gene.activities)
so <- NormalizeData(
  object = so,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(so$nCount_RNA)
)

DefaultAssay(so) <- "RNA"
FeaturePlot(so, features = c("CCR7", "SELL", "IL7R", "S1PR5", "ZEB2", "ID3", "TCF7", "CD8A", "CD8B", "CD4"), max.cutoff = "q95", 
            reduction = "peaks.umap") 
  




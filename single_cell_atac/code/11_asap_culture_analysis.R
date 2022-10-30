library(Seurat)
library(Signac)
library(BuenColors)
library(dplyr)
library(data.table)
library(Matrix)
library(viridis)

# Import and set the table
sc <- fread("../data/asap_culture_data/TEM_culture-ATAC_hg38_v20-mtMask.singlecell.csv")
mito <- readRDS("../data/asap_culture_data/TEM_culture-ATAC_hg38_v20-7076.rds")
h5 <- Read10X_h5("../data/asap_culture_data/TEM_culture-ATAC_hg38_v20.peakbycell.h5")

mito <- mito %>% dplyr::filter(cov7076 >= 10 & cov7076 < 200)
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
mdf <- mdf %>% dplyr::filter(cell %in% colnames(tags))



# Do the Seurat object now
chrom_assay <- CreateChromatinAssay(
  counts = h5[grepl("^chr", rownames(h5)),mdf$cell],
  sep = c(":", "-"),
  #genome = 'hg38',
  fragments = '../../../m7076-large-data-files/fragments-files/TEM_culture-ATAC_hg38_v20.fragments.tsv.gz',
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

so$frip <- so$peak_region_fragments/so$passed_filters
so$tsspct <- so$TSS_fragments/so$passed_filters

qplot(log10(so@meta.data$passed_filters), so$frip)
#so <- subset(so, subset = passed_filters > 1000 & so$frip > 0.4)

# Add ADT data
so[["ADT"]] <- CreateAssayObject(counts = tags[,colnames(so)])
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
  RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:15, reduction.name = "peaks.umap")
so <- so %>% FindNeighbors( reduction = 'lsi', dims = 2:15, graph.name = 'peaks_lsi')
so <- so %>%  FindClusters( graph.name = "peaks_lsi", resolution = 0.7)

DimPlot(object = so, label = TRUE, reduction = "peaks.umap") 

sc <- so$seurat_clusters
so$caleb_cluster <- case_when(
  sc == "0" ~ "CD8 MPEClike",
  sc == "1" ~ "CD8 SLEClike",
  sc == "2" ~ "CD4 Cytotoxic",
  sc == "3" ~ "CD4 Effector",
  sc == "4" ~ "CD8 Terminal Effector",
  sc == "5" ~ "CD8 MPEClike",
  sc == "6" ~ "CD4 Effector",
  sc == "7" ~ "CD8 SLEClike"
)

clusters_umap <- DimPlot(object = so, label = FALSE, reduction = "peaks.umap", group.by = "caleb_cluster",
                         pt.size = 0.4) +
  scale_color_manual(values = jdb_palette("brewer_spectra")[c(1,2,7,9,8)]) +
  theme_void() + theme(legend.position = "none") +
  ggtitle("")
cowplot::ggsave2(clusters_umap, file = paste0("../output/cultureplots/clusters_umap_viz.png"), width = 4, height = 4, dpi = 500)


# See enrichment of A/G allele in clusters
so@meta.data %>%
  group_by(caleb_cluster) %>%
  mutate(hi7076 = af7076 > 0.9, low7076 = af7076 < 0.1) %>%
  summarize(
    ratio = sum(low7076)/sum(hi7076),
    n_hi = sum(hi7076),
    n_low = sum(low7076),
    proportion = sum(hi7076) / (sum(hi7076) + sum(low7076)),
    count = sum(hi7076) +  sum(low7076),
    mean_heteroplasmy = round(mean(af7076),3)) %>%
  arrange(desc(ratio)) %>% mutate(rank = 1:n()) -> outdf

rest_n <- sum(outdf$count[2:5])
rest_prop <- sum(outdf$proportion[2:5]*outdf$count[2:5])/rest_n

binom.test(p = rest_prop, x = outdf$n_hi[1], n = outdf$count[1])

outdf$caleb_cluster_order <- factor(outdf$caleb_cluster, levels = outdf$caleb_cluster)
p2 <- ggplot(outdf, aes(x = caleb_cluster_order, y = ratio, color = caleb_cluster)) +
  geom_point(size = 3) +
  scale_color_manual(values = jdb_palette("brewer_spectra")[c(1,2,7,9,8)]) +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.9, 1.35))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(x = "T cell cluster", y = "")

p2
cowplot::ggsave2(p2, file = "../output/cultureplots/rank_order_plot_ratio.pdf", width = 1.4, height = 1.6)


DefaultAssay(so) <- "ADT"
FindAllMarkers(so, only.pos = TRUE)

library(viridis)
so$wildtype_weight <- Nebulosa:::calculate_density(1-(so$af7076 ), (so@reductions$peaks.umap@cell.embeddings[,c(1,2)]),
                                                   method = "wkde", adjust = 2)
so$klrg1_weight <- Nebulosa:::calculate_density(so@assays$ADT@data["KLRG1",], (so@reductions$peaks.umap@cell.embeddings[,c(1,2)]),
                                                method = "wkde", adjust = 2)
so$il7r_weight <- Nebulosa:::calculate_density(so@assays$ADT@data["IL7R",], (so@reductions$peaks.umap@cell.embeddings[,c(1,2)]),
                                               method = "wkde", adjust = 2)

# Look at heteroplasmy
plot_het_umap_smooth <- ggplot(cbind(so@meta.data,
                                     so@reductions$peaks.umap@cell.embeddings[,c(1,2)])%>% arrange(wildtype_weight),
                               aes(x = UMAP_1, y = UMAP_2, color = wildtype_weight)) +
  geom_point(size = 1) +
  scale_color_gradientn(colors = jdb_palette("flame_flame")) + 
  theme_void() +
  theme(legend.position = "none")
cowplot::ggsave2(plot_het_umap_smooth, file = paste0("../output/cultureplots/smoothed_umap_viz.png"), width = 4, height = 4, dpi = 500)

# Plot markers nebulosa
plot_klrg1_smooth <- ggplot(cbind(so@meta.data,
                                     so@reductions$peaks.umap@cell.embeddings[,c(1,2)])%>% arrange(klrg1_weight),
                               aes(x = UMAP_1, y = UMAP_2, color = klrg1_weight)) +
  geom_point(size = 1) +
  scale_color_viridis() +
  theme_void() +
  theme(legend.position = "none")
cowplot::ggsave2(plot_klrg1_smooth, file = paste0("../output/cultureplots/smoothedKLRG1_umap_viz.png"), width = 4, height = 4, dpi = 500)

plot_il7r_smooth <- ggplot(cbind(so@meta.data,
                                     so@reductions$peaks.umap@cell.embeddings[,c(1,2)])%>% arrange(il7r_weight),
                               aes(x = UMAP_1, y = UMAP_2, color = il7r_weight)) +
  geom_point(size = 1) +
  scale_color_viridis() +
  theme_void() +
  theme(legend.position = "none")
cowplot::ggsave2(plot_il7r_smooth, file = paste0("../output/cultureplots/smoothedIL7R_umap_viz.png"), width = 4, height = 4, dpi = 500)


# CD4/CD8 ratio
so$CD4CD8ratio <- log2((so@assays$ADT@counts["CD8A",]+1)/(so@assays$ADT@counts["CD4",]+1))
so$smooth_cd4cd8 <- Nebulosa:::calculate_density(so$CD4CD8ratio, (so@reductions$peaks.umap@cell.embeddings[,c(1,2)]),
                                               method = "wkde", adjust = 0.5)
umapCD4CD8 <- ggplot(cbind(so@meta.data,
                                     so@reductions$peaks.umap@cell.embeddings[,c(1,2)])%>% arrange((abs(CD4CD8ratio))),
                               aes(x = UMAP_1, y = UMAP_2, color = CD4CD8ratio)) +
  geom_point(size = 1) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) + 
  theme_void() +
  theme(legend.position = "none")
cowplot::ggsave2(umapCD4CD8, file = paste0("../output/cultureplots/cd4cd8_umap_viz.png"), width = 4, height = 4, dpi = 500)


# regulat hetero
plot_het_umap_reg <- ggplot(cbind(so@meta.data,
                                  so@reductions$peaks.umap@cell.embeddings[,c(1,2)])%>% arrange(desc(af7076)),
                            aes(x = UMAP_1, y = UMAP_2, color = af7076)) +
  geom_point(size = 1) +
  scale_color_viridis()+
  theme_void() +
  theme(legend.position = "none")
cowplot::ggsave2(plot_het_umap_reg, file = paste0("../output/cultureplots/regularhet_umap_viz.png"), width = 4, height = 4, dpi = 500)

###
DefaultAssay(so) <- "ADT"
FeaturePlot(so, features = c("CD4CD8ratio", "klrg1_weight", "il7r_weight",  "wildtype_weight", "ITGB1"),
            max.cutoff = "q99", min.cutoff = 'q05',
            reduction = "peaks.umap") &
  scale_color_viridis()


sso2 <- subset(so, subset = seurat_clusters %in% c(0,1,4,5,7)) # c(2,3,6)
cordf <- data.frame(
  cor = cor(t(data.matrix(so2@assays$ADT@data)), so2$af7076))
cordf$p_value <- sapply(1:dim(so2@assays$ADT@data)[1], function(i){
  cor.test(so2@assays$ADT@data[i,], so2$af7076)$p.value
})
cordf <- cordf %>% arrange(desc(cor))
cordf$p_adj <- p.adjust(cordf$p_value)
cordf

summary(so@meta.data$af7076)


if(FALSE){
  
  # Now create the Seurat object
  library(EnsDb.Hsapiens.v86)
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  genome(annotations) <- "hg38"
  seqlevelsStyle(annotations) <- 'UCSC'
  
  # Annotate
  annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
  Annotation(so) <- annotations
  
  # View gene scores
  DefaultAssay(so) <- "peaks"
  gene.activities <- GeneActivity(so)
  so[['RNA']] <- CreateAssayObject(counts = gene.activities)
  so <- NormalizeData(
    object = so,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(so$nCount_RNA)
  )
  
  DefaultAssay(so) <- "RNA"
  fma <- FindAllMarkers(so, only.pos = TRUE)
  FeaturePlot(so, features = c("CCR7", "SELL", "IL7R", "IKZF2", "ZEB2", "IFNG", "TCF7", "CD8A", "CD8B", "CD4"), max.cutoff = "q95", 
              reduction = "peaks.umap", order = TRUE) 
  
  DefaultAssay(so) <- "RNA"
  plot_one_gs <- function(gene){
    p1 <- FeaturePlot(so, features = gene, reduction = "peaks.umap", pt.size = 0.75,
                min.cutoff = "q1", max.cutoff = "q95", order = TRUE) +
      theme_void() + theme(legend.position = "none") + ggtitle("")
    cowplot::ggsave2(p1, file = paste0("../output/cultureplots/", gene, "_viz.png"), width = 4, height = 4, dpi = 500)
  }
  plot_one_gs("IKZF2")
  plot_one_gs("ZEB2")
  plot_one_gs("IFNG")
  plot_one_gs("KLRD1")
  
}




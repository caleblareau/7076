library(data.table)
library(dplyr)
library(BuenColors)
"%ni%" <- Negate("%in%")

dates <- c(20201211, 20210115, 20210218, 20210329, 20210511)

date0 <- as.numeric(as.Date(as.character(20201211), format = "%Y%m%d"))
import_day <- function(date){
  
  # Format date
  Day_n <- as.numeric(as.Date(as.character(date), format = "%Y%m%d")) - date0
  
  # Import data
  multiplet <- fread(paste0("../data/atac_data/ATS_",date,"_C1_MultipletBarcodes.txt"), header = FALSE)[[1]]
  m7076 <- readRDS(paste0("../data/mito_data/m7076/ATS_PBMC_",date,"_C1_hg38_v20-7076.rds"))
  mdf <- merge(fread(paste0("../data/atac_data/ATS_",date,"_C1_refmapped.csv.gz")), m7076, by.x = "cb", by.y = "cell")
  
  # Now do some dplyr stuff
  mdf %>% filter(cb %ni% multiplet & cov7076 > 10) %>%
    mutate(Day = Day_n)
  
  mdf$multiplet <- mdf$cb %in% multiplet
  
  mdf %>% mutate(Day = Day_n)
}

# Do some QC
lapply(dates, import_day) %>% rbindlist() %>% data.frame() -> full_df
table(full_df$Day)
full_df <- full_df %>% filter(!multiplet)
dim(full_df)

# Create an ideal annotation
full_df$annotation <- case_when(
  full_df$predicted.celltype.l1 == "Mono" ~ "Mono",
  full_df$predicted.celltype.l1 == "B" ~ "Bcell",
  full_df$predicted.celltype.l2 == "CD4 TCM" ~ "CD4tcm",
  full_df$predicted.celltype.l2 == "CD8 TEM" ~ "CD8tem",
  full_df$predicted.celltype.l2 == "CD8 Naive" ~ "CD8naive",
  full_df$predicted.celltype.l2 == "CD4 Naive" ~ "CD4naive",
  full_df$predicted.celltype.l1 == "NK" ~ "NK",
  TRUE ~ "other"
)
table(full_df$annotation)


so <- readRDS("../../../7076-large-data-files/seurat-projections/ATS_20201211_C1_refumap.rds")

df <- full_df %>% filter(orig.ident == "donor3_ATS_20201211_C1")
rownames(df) <- df$cb
so <- AddMetaData(so, df)

so$X7076G <- so$af7076 > 0.9
so2 <- subset(so, subset = (af7076 > 0.9 |af7076 < 0.1) & predicted.celltype.l2 == "CD4 Naive")
so.tem <- subset(so, subset = (predicted.celltype.l2 == "CD4 Naive"))

library(Seurat)
so2 <- ScaleData(so2, assay = "ACTIVITY")
diff_gene_scores <- FindMarkers(so2, assay = "ACTIVITY", group.by = "X7076G",
                                ident.1 = "TRUE", ident.2 = "FALSE", logfc.threshold  = 0.1)
diff_gene_scores

cordf <- data.frame(
  cor = cor(t(data.matrix(so.tem@assays$ACTIVITY@data)), so.tem$af7076)
) %>% arrange(desc(cor))

head(cordf)

cordf[!is.na(cordf$cor),, drop = FALSE] %>% tail(20)

library(data.table)
library(dplyr)
library(BuenColors)
library(viridis)
"%ni%" <- Negate("%in%")

theme_blank <- function(...) {
  ret <- theme_bw(...)
  ret$line <- element_blank()
  ret$rect <- element_blank()
  ret$strip.text <- element_blank()
  ret$axis.text <- element_blank()
  ret$plot.title <- element_blank()
  ret$axis.title <- element_blank()
  ret$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")
  ret
}

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

full_df %>%
  mutate(hi7076 = af7076 > 0.9, low7076 = af7076 < 0.1) %>%
  group_by(predicted.celltype.l2) %>%
  summarize(n = n(), hi7076 = sum(hi7076), low7076 = sum(low7076)) %>%
  filter(n > 50) %>%
  mutate(ratio = hi7076/low7076) %>%
  arrange(desc(ratio)) %>%
  ungroup() %>%
  mutate(rank = 1:n()) -> out_df

out_df$fill_color <- case_when(
  out_df$predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono") ~ "Mono",
  out_df$predicted.celltype.l2 %in% c("B memory", "B naive", "B intermediate")~ "Bcell",
  out_df$predicted.celltype.l2 == "CD4 TCM" ~ "CD4tcm",
  out_df$predicted.celltype.l2 == "CD8 TEM" ~ "CD8tem",
  out_df$predicted.celltype.l2 == "CD8 Naive" ~ "CD8naive",
  out_df$predicted.celltype.l2 == "CD4 Naive" ~ "CD4naive",
  out_df$predicted.celltype.l2 %in% c("NK", "NK_CD56bright") ~ "NK",
  TRUE ~ "other"
)

p2 <- ggplot(out_df, aes(x = rank, y = ratio, fill = fill_color)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = jdb_palette("corona")[c(5,2:4,1,6,7,8)]) +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")

overall_mean <- sum(out_df$hi7076)/(sum(out_df$hi7076) + sum(out_df$low7076))
out_df$p_value <- pbinom(out_df$hi7076, out_df$hi7076 + out_df$low7076,  prob = overall_mean, lower.tail = TRUE)

ggplot(out_df, aes(x = ratio, y = -log10(p_value), color = fill_color)) +
  geom_point() +
  scale_color_manual(values = jdb_palette("corona")[c(5,2:4,1,6,7,8)]) +
  pretty_plot(fontsize = 7) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none")

cowplot::ggsave2(p2, file = "../output/rank_order_plot-recip.pdf", width = 1.6, height = 1)

# Try making a UMAP viz first
library(viridis)

ggplot(shuf(full_df), aes(x = refUMAP_1, y = refUMAP_2, color = af7076)) +
  scale_color_viridis() + theme_blank() +
  theme(legend.position = "none") +
  geom_density_2d() +
  facet_wrap(~af7076 > 0.95 )

plot_het_umap <- ggplot(shuf(full_df), aes(x = refUMAP_1, y = refUMAP_2, color = af7076)) +
  geom_point(size = 0.1) + scale_color_viridis() + theme_blank() +
  theme(legend.position = "none")

cowplot::ggsave2(plot_het_umap, width = 5, height = 5, filename = "../output/heteroplasmy_umap.png")

full_df$af_weight <- Nebulosa:::calculate_density(1-(full_df$af7076 ), cbind(full_df$refUMAP_1, full_df$refUMAP_2),
                                                  method = "wkde", adjust = 5)
plot_het_umap_smooth <- ggplot(full_df %>% arrange(af_weight),
                               aes(x = refUMAP_1, y = refUMAP_2, color = af_weight)) +
  geom_point(size = 0.1) +
  scale_color_gradientn(colors = jdb_palette("flame_flame")) + theme_blank() +
  theme(legend.position = "none")

cowplot::ggsave2(plot_het_umap_smooth, width = 5, height = 5, filename = "../output/heteroplasmy_smooth_umap.png")

full_df$high <- full_df$af7076 > 0.95
full_df$low <- full_df$af7076 < 0.05

plot_ct_umap <- ggplot(full_df , aes(x = refUMAP_1, y = refUMAP_2, color = annotation)) + 
  geom_point(size = 0.1) +
  scale_color_manual(values = jdb_palette("corona")[c(5,2:4,1,6,7,8)]) +
  theme_blank() + theme(legend.position = "none")
cowplot::ggsave2(plot_ct_umap, width = 5, height = 5, filename = "../output/celltype_umap.png")



ggplot(full_df, aes(x = refUMAP_1, y = refUMAP_2, color = predicted.celltype.l2)) + 
  geom_point() + scale_color_manual(values = jdb_palette("corona")) +
  facet_wrap(~Day)

# Look at overall abundance
p1 <- ggplot(full_df, aes(x = af7076*100)) +
  geom_histogram(bins = 11, fill = "lightgrey", color = "black") +
  pretty_plot(fontsize = 7) + scale_y_continuous(expand = c(0,0)) +
  L_border() + labs(x = "m.7076A>G heteroplasmy (%)", y = "single cell count")
p1
cowplot::ggsave2(p1, file = "../output/basic_histo_panelB.pdf", width = 1.5, height = 1.1)


pSupplement <- ggplot(full_df %>%
                        group_by(predicted.celltype.l2) %>% mutate(count = n()) %>% filter(count > 100), aes(x = af7076*100)) +
  geom_histogram(bins = 11, fill = "lightgrey", color = "black") +
  pretty_plot(fontsize = 7) + scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~predicted.celltype.l2, scale = "free_y") +
  labs(x = "m.7076A>G heteroplasmy (%)", y = "single cell count")

cowplot::ggsave2(pSupplement, file = "../output/supplement_histogram.pdf", 
                 width = 4, height = 3)

ggplot(full_df  %>% filter(annotation != "other"), aes(x = af7076*100)) +
  geom_histogram(bins = 11, fill = "lightgrey", color = "black") +
  pretty_plot(fontsize = 7) + scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~annotation, scale = "free_y") +
  labs(x = "m.7076A>G heteroplasmy (%)", y = "single cell count")

full_df  %>% filter(annotation != "other") %>%
  group_by(Day, annotation)  %>%
  mutate(af7076 = af7076*100) %>%
  summarize(maf = mean(af7076), sem = sd(af7076)/sqrt(n())) -> long_data_df

long_data_df %>% 
  ggplot(aes(x = Day, y = maf, color = annotation)) +
  geom_point(size = 0.3) + geom_line() +
  geom_errorbar(aes(ymin=maf-sem, ymax=maf+sem), width=30,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 47.3, linetype = 2, color = "black") + 
  pretty_plot(fontsize = 6) + labs(x = "Day after initial sample", y = "Mean m.7076A>G heteroplasmy (%)") +
  theme(legend.position = "bottom") + L_border() +
  scale_color_manual(values = jdb_palette("corona")[c(5,2:4,1,6,7)]) -> plot_long
cowplot::ggsave2(plot_long, width = 1.5, height = 2, filename = "../output/longitudinal_het.pdf")

long_data_df %>% 
  ggplot(aes(x = Day, y = maf, color = annotation)) +
  geom_point(size = 0.3) + geom_line() +
  geom_errorbar(aes(ymin=maf-sem, ymax=maf+sem), width=30,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 47.3, linetype = 2, color = "black") + 
  pretty_plot(fontsize = 6) + labs(x = "Day after initial sample", y = "Mean m.7076A>G heteroplasmy (%)") +
  theme(legend.position = "bottom") + L_border() + scale_y_continuous(limits = c(0,55))+
  scale_color_manual(values = jdb_palette("corona")[c(5,2:4,1,6,7)]) -> plot_long
cowplot::ggsave2(plot_long, width = 1.5, height = 2, filename = "../output/longitudinal_het_from0.pdf")


# For supplement
smoothScatter(full_df$cov7076,full_df$af7076*100, 
              colramp = colorRampPalette(c("white", jdb_palette("solar_rojos"))), pch = NA,
              xlab = "Coverage at m7076A>G", ylab = "Single-cell heteroplasmy") 

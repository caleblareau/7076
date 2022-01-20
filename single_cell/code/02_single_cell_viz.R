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
  
}

lapply(dates, import_day) %>% rbindlist() %>% data.frame() -> full_df
table(full_df$Day)

p1 <- ggplot(full_df, aes(x = af7076*100)) +
  geom_histogram(bins = 11, fill = "lightgrey", color = "black") +
  pretty_plot(fontsize = 7) + scale_y_continuous(expand = c(0,0)) +
  L_border() + labs(x = "m.7076A>G heteroplasmy (%)", y = "single cell count")
cowplot::ggsave2(p1, file = "../output/basic_histo_panelB.pdf", width = 1.5, height = 1.1)

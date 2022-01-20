library(data.table)
library(dplyr)
library(BuenColors)
files <- list.files("../data/mito_data/var_stats/", full.names = TRUE)
lapply(files, function(f){
  fread(f) %>%
    filter(n_cells_conf_detected >= 1 & strand_correlation > 0.65 & vmr > 0.01 & mean_coverage > 20) %>%
    arrange(desc(mean))
}) -> list_of_var_stats

called_variants <- Reduce(intersect, lapply(list_of_var_stats, function(x) x$variant))
var_df <- fread(files[1]) %>%
  data.frame() %>% filter(variant %in% called_variants) %>%
  arrange(desc(mean)) %>% 
  mutate(rank = 1:n())
dim(var_df)
p1 <- ggplot(var_df, aes(x = rank, y = mean*100)) +
  geom_point(size = 0.3) + 
  pretty_plot(fontsize = 7) + labs(x = "Rank ordered mutations (n = 183)", y = "% pseduobulk heteroplasmy ") +
  L_border()
median(var_df$mean)*100
cowplot::ggsave2(p1, file = "../output/rank_order_muts.pdf", width = 1.5, height = 1.5)

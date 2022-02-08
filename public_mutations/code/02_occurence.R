library(dplyr)
library(BuenColors)
library(ggforce)
library(data.table)

source("01_functions.R")

# Set up coordinates
df <- data.frame(pos = 1:16569, angle = seq(0, (2*pi), by = (2*pi)/(16569-1)))

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("m", annotations$Position, annotations$Reference, ">", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

# Look at specific amino acid changes
annotations %>% filter(syn_annotation == "WCF_to_WCF") %>% pull(Amino.acids) %>% table()
annotations %>% filter(syn_annotation == "Wobble_to_Wobble") %>% pull(Amino.acids) %>% table()
annotations %>% filter(syn_annotation == "WCF_to_Wobble") %>% pull(Amino.acids) %>% table()
annotations %>% filter(syn_annotation == "Wobble_to_WCF") %>% pull(Amino.acids) %>% table()

# Sanity check
annotations %>% filter(Position == 7076) 

pos_count_df <- annotations %>% group_by(Position, syn_annotation) %>% summarize(count = n()) %>%
  filter(syn_annotation != "other")
mmdf <- merge(df, pos_count_df, by.x = "pos", by.y = 'Position')
mmdf$syn_annotation_int <- as.integer(as.factor(mmdf$syn_annotation))
count_me(annotations)

library(ggbeeswarm)
p1 <- ggplot(mmdf, aes(x = syn_annotation_int, y = pos, color = syn_annotation))+
  geom_quasirandom(size = 0.1) + coord_polar(theta = "y") +
  theme_void() + 
  xlim(-5,5) + ylim(1, 16569) +
  scale_color_manual(values = jdb_palette("corona"))

cowplot::ggsave2(p1, file = "../output/possible_mutations.pdf", width = 3, height = 3)

library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)

source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("m", annotations$Position, annotations$Reference, ">", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

# Import Phewas associations
assoc <- readxl::read_xlsx("../data/austim_Wang2022NatureComms_41467_2022_30805_MOESM3_ESM.xlsx", 1) %>% data.frame()

assoc$mutation <- paste0("m", assoc$mtDNA.position, toupper(assoc$Reference.allele),">",toupper(assoc$Variant.allele))

mdf <- merge(assoc, annotations, by = "mutation") 


# Count and compare
all_df_count <- rbind(
  count_me(mdf  ) %>% mutate(what = "autism"),
  count_me(annotations) %>% mutate(what = "expected")
)

sum(all_df_count$count[1:4])

library(ggbeeswarm)
p1bar <- ggplot(all_df_count %>% filter(what == "autism"), aes(x = syn_annotation, y = perc, fill = syn_annotation)) +
  geom_bar( stat = "identity", position = "dodge", color = "black",) +
  scale_fill_manual(values = jdb_palette("corona")[1:4] )+
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), limits = c(0,45))+ 
  labs(x = "", y = "% of synonymous mtDNA variants") + theme(legend.position = "none")
p1bar
cowplot::ggsave2(p1bar, file = "../output/freqs_syn_type_autsm.pdf", width = 2, height = 1.5)

c(prop.test(all_df_count$count[1], sum(all_df_count$count[1:4]), p = all_df_count$perc[5] / 100)$p.value,
  prop.test(all_df_count$count[2], sum(all_df_count$count[1:4]), p = all_df_count$perc[6] / 100)$p.value,
  prop.test(all_df_count$count[3], sum(all_df_count$count[1:4]), p = all_df_count$perc[7] / 100)$p.value,
  prop.test(all_df_count$count[4], sum(all_df_count$count[1:4]), p = all_df_count$perc[8] / 100)$p.value)

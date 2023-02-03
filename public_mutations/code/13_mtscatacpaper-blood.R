library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
library(stringr)
source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0( annotations$Position, annotations$Reference, ">", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

# Import Phewas associations
vars <- rownames(readRDS("../data/filteredCD34_mgatk_calls.rds"))
rare_vars <- rownames(readRDS("../data/mtscatac_fig6_rare_cd34_vars.rds"))
mdf <- merge(data.frame(mutation = c(vars, rare_vars)), annotations, by = "mutation") 

process_stats_syn_type <- function(df){
  cm_a <- count_me(annotations)
  cm_t <- count_me(df)
  cm_t$all_perc <- cm_a$perc
  cm_t$pvalue <- c(prop.test(cm_t$count[1], sum(cm_t$count), p = cm_a$perc[1] / 100)$p.value,
                   prop.test(cm_t$count[2], sum(cm_t$count), p = cm_a$perc[2] / 100)$p.value,
                   prop.test(cm_t$count[3], sum(cm_t$count), p = cm_a$perc[3] / 100)$p.value,
                   prop.test(cm_t$count[4], sum(cm_t$count), p = cm_a$perc[4] / 100)$p.value)
  cm_t
  
}
process_stats_syn_type(annotations) 
process_stats_syn_type(mdf ) 

# Count and compare
all_df_count <- rbind(
  count_me((mdf %>% filter(donor_type == "s1") )  ) %>% mutate(what = "autism"),
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

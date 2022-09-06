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

# Merge the two and look at FDR < 0.05 hits (how significance has been previously defined)
mdf <- merge(assoc, annotations, by = "mutation") 
mdf %>%
  mutate(fdr_sig = FDR < 0.05) %>%
  group_by(mutation, syn_annotation, Consequence) %>% summarize(count = sum(fdr_sig)) %>%
  mutate(is_sig = count > 0) %>%
  group_by(is_sig, syn_annotation, Consequence) %>% summarize(count = n())



# Count and compare
all_df_count <- rbind(
  count_me(mdf  ) %>% mutate(what = "expected"),
  count_me(annotations) %>% mutate(what = "autism")
)

library(ggbeeswarm)
p1bar <- ggplot(all_df_count, aes(x = syn_annotation, y = perc, fill = syn_annotation)) +
  geom_bar(data = cm_c, stat = "identity", position = "dodge", color = "black",) +
  scale_shape_manual(values=c(0, 1, 2,3))+
  scale_fill_manual(values = jdb_palette("corona")[1:4] )+
  pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(expand = c(0,0), limits = c(0,45))+ 
  labs(x = "", y = "% of synonymous mtDNA variants") + theme(legend.position = "none")
p1bar
cowplot::ggsave2(p1bar, file = "../output/freqs_syn_type.pdf", width = 2, height = 1.5)


kdf20 <- kdf %>%
  filter(FDR < 0.2) 


# Look at what recurrent phenotypes we have
kdf05[,c("mutation", "Consequence")]%>% group_by(mutation, Consequence) %>%
  summarize(count =n ()) %>% pull(Consequence) %>% table()

# Look at the patterns of mutations of interest
kdf20[grepl("\\*",kdf20$Reference.tRNA) & grepl("near",kdf20$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[grepl("\\*",kdf05$Reference.tRNA) & grepl("near",kdf05$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[grepl("\\*",kdf05$Variant.tRNA) & grepl("near",kdf05$Reference.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[grepl("near",kdf05$Reference.tRNA) & grepl("near",kdf05$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)


kdf05[grepl("\\*",kdf05$Reference.tRNA) & grepl("\\*",kdf05$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[,c("Trait", "Consequence")]%>% filter(Consequence == "synonymous_variant") %>%
  group_by(Trait) %>% summarize(count = n()) %>% arrange(desc(count))


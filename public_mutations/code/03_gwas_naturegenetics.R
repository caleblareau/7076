library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)

source("01_functions.R")

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("m", annotations$Position, annotations$Reference, ">", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

annotations[,c("mutation","Position", "Reference", "Variant", "Symbol", "Codons", "Reference.tRNA","Variant.tRNA","syn_annotation")] %>%
  filter(!(syn_annotation == "other")) %>%
  write.table("../output/syn-table-for-ryan.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Import Phewas associations
assoc <- readxl::read_xlsx("../data/MitoPhewas_associations.xlsx", 1) %>% data.frame()
assoc$mutation <- paste0("m", assoc$Position, toupper(assoc$Allele1),">",toupper(assoc$Allele2))
assoc_2 <- assoc
assoc_2 <- assoc_2 %>% filter(mutation %in% (annotations %>%pull(mutation)))
assoc_2$is_syn <- assoc_2$mutation %in% (annotations %>% filter(syn_annotation != "other") %>%
                                           pull(mutation))
assoc_2 %>%
  mutate(fdr_sig = FDR < 0.05) %>%
  group_by(mutation, is_syn) %>% summarize(count = sum(fdr_sig)) %>%
  mutate(is_sig = count > 0) %>%
  group_by(is_sig, is_syn) %>% summarize(count = n())

assoc_2 %>%
  mutate(fdr_sig = FDR < 0.05) %>%
  group_by(fdr_sig, is_syn) %>% summarize(count = sum(fdr_sig)) 


# Merge the two and look at FDR < 0.05 hits (how significance has been previously defined)
mdf <- merge(assoc, annotations, by = "mutation") %>% arrange(FDR)
mdf %>%
  filter( FDR < 0.05) %>%
  group_by( syn_annotation) %>%
  summarize(count = n())


  mutate(is_sig = count > 0) %>%
  group_by(is_sig, syn_annotation, Consequence) %>% summarize(count = n())



prop.test(4,170, 0.0421)
# Count and compare
count_me(mdf %>% filter(FDR < 0.05)  )
count_me(annotations)
count_me(annotations %>% filter(Disease != "-"))

mdf %>% filter(FDR < 0.05) %>% filter(syn_annotation!= "other") %>%
  arrange(desc(syn_annotation)) %>% filter(abs(Effect.size)>1) %>% 
  arrange(desc(Effect.size))
  pull(Effect.size) %>% summary()

annotations %>% filter(Disease != "-") %>% filter(syn_annotation != "other") %>%
  arrange(syn_annotation)

kdf <- mdf[,c("Position.x","mutation","Symbol","Trait", "FDR", "Codons","Consequence", "Reference.tRNA", "Variant.tRNA", "Effect.size", "syn_annotation")] 
kdf05 <- kdf %>%
  filter(FDR < 0.05) 


kdf05 %>% group_by(Symbol, mutation) %>% summarize(count = n()) %>%
  group_by(Symbol) %>% summarize(count = n()) %>%
  ungroup() %>% mutate(perc = count/sum(count)*100)

# Make a circos plot
gwas_plot <- kdf05 %>% group_by(Consequence, Symbol, Position.x, mutation, syn_annotation) %>% summarize(count = n()) %>%
  filter(syn_annotation != "other")
gwas_plot$syn_annotation_int <- as.integer(as.factor(gwas_plot$syn_annotation))

library(ggbeeswarm)
p1 <- ggplot(gwas_plot, aes(x = syn_annotation_int, y = Position.x, color = syn_annotation))+
  geom_quasirandom(size = 0.1) + coord_polar(theta = "y") +
  theme_void() + 
  xlim(-5,5) + ylim(1, 16569) +
  scale_color_manual(values = jdb_palette("corona"))

#cowplot::ggsave2(p1, file = "../output/gwas_mutations.pdf", width = 3, height = 3)




kdf20 <- kdf %>%
  filter(FDR < 0.2) 


# Look at what recurrent phenotypes we have
kdf05[,c("mutation", "Consequence")]%>% group_by(mutation, Consequence) %>%
  summarize(count =n ()) %>% pull(Consequence) %>% table()

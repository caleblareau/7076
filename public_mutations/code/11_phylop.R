library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
library(stringr)
source("01_functions.R")

 # Import annotation
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("m", annotations$Position, annotations$Reference, ">", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

####

annotations %>% filter(syn_annotation != "other") %>%
  group_by(Position, syn_annotation) %>% summarize(count = n()) %>%
  filter(count == 3) %>% pull(syn_annotation) %>% table()
  
phylop20 <- fread("../data/phylop20.chrM.txt")
vec <- phylop20[[4]]; names(vec) <- as.character(phylop20[[3]])
annotations$phylop <- vec[as.character(annotations$Position)]

annotations %>% filter(syn_annotation != "other") %>%
  group_by(Position, syn_annotation) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(Position) %>% mutate(total = sum(count)) %>%
  mutate(all_synon = total == 3) %>% group_by(Position, all_synon) %>% 
  summarize(class = paste0(syn_annotation, collapse = "-")) -> anno_me
anno_me$class2 <- case_when(
  anno_me$all_synon ~ anno_me$class,
  !anno_me$all_synon ~ "nonAllSynon"
)
anno_me$phylop <- vec[as.character(anno_me$Position)]
anno_me %>%
  ggplot(aes(x = class2, y = phylop)) + pretty_plot() +
  geom_boxplot(outlier.shape = NA, width = 0.5) + labs(x = "Codon class") +
  pretty_plot(fontsize = 7) + L_border() + coord_flip(ylim = c(-11, 2)) -> bah
cowplot::ggsave2(bah, file = "../output/boxplot.pdf", width = 3.1, height = 1.8)

wilcox.test(
  anno_me$phylop[anno_me$class == "WCF_to_Wobble"],
  anno_me$phylop[anno_me$class == "Wobble_to_WCF-Wobble_to_Wobble"]
) %>% str()

wilcox.test(
  anno_me$phylop[anno_me$class2 == "WCF_to_Wobble"],
  anno_me$phylop[anno_me$class2 == "nonAllSynon"]
) %>% str()


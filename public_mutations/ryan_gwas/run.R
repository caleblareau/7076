library(data.table)
library(dplyr)
library(BuenColors)
library(stringr)

quant <- fread("mt_syn_quant.tsv")
quant$position <- as.numeric(as.character(str_split_fixed(quant$geno, "-", 4)[,2]))
p1 <- ggplot(quant , aes(x = position, y = -log10(p + 1e-196), color = syn_annotation)) +
  geom_point(size = 0.5) + facet_grid( ~ syn_annotation) +
  pretty_plot(fontsize = 7) + scale_color_manual(values = jdb_palette("corona")) +
  theme(legend.position = "none") + labs(x = "Position in mtDNA genome", y = "-log10 p-value")
cowplot::ggsave2(p1, file = "association_binary.pdf", width = 4.2, height = 1.3)


quant %>% arrange(pValue) %>% filter(homoplasmic_MAF < 0.05) %>% 
  head(20)
quant %>% filter(p < 1e-50) %>% arrange(syn_annotation)
quant %>% filter(pheno == "Reticulocyte count") %>%
  mutate(up = beta > 0) %>% 
  group_by(up, syn_annotation) %>% summarize(count = n())

#########3


bin <- fread("mt_syn_binary.tsv") %>% filter(model == "Ryan")
bin$position <- as.numeric(as.character(str_split_fixed(bin$genotype, "-", 4)[,2]))
p2 <- ggplot(bin , aes(x = position, y = -log10(pValue + 1e-196), color = syn_annotation)) +
  geom_point(size = 0.5) + facet_grid( ~ syn_annotation) +
  pretty_plot(fontsize = 7) + scale_color_manual(values = jdb_palette("corona")) +
  theme(legend.position = "none") + labs(x = "Position in mtDNA genome", y = "-log10 p-value")
cowplot::ggsave2(p2, file = "association_binary.pdf", width = 4.2, height = 1.3)




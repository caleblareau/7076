library(data.table)
library(dplyr)
library(BuenColors)
library(stringr)

quant <- fread("MT_quant_SYN_p0001.tsv") 
quant$position <- as.numeric(as.character(str_split_fixed(quant$geno, "-", 4)[,2]))
p1 <- ggplot(quant  , aes(x = position, y = -log10(p + 1e-196), color = syn_annotation)) +
  geom_point(size = 0.5) + facet_grid( ~ syn_annotation) +
  pretty_plot(fontsize = 7) + scale_color_manual(values = jdb_palette("corona")) +
  theme(legend.position = "none") + labs(x = "Position in mtDNA genome", y = "-log10 p-value")
cowplot::ggsave2(p1, file = "plots/association_quant_maf05-updated.pdf", width = 4, height = 1.2)


#########3


bin <- fread("MT_binary_SYN_p0001.tsv") %>% filter(model == "Ryan2")
bin$position <- as.numeric(as.character(str_split_fixed(bin$genotype, "-", 4)[,2]))
bin[-25,c("genotype", "pValue", "TT", "Symbol", "syn_annotation", "phenotype")]  %>% filter( pValue < 10^-4) %>% data.frame() %>% head(40)
bin[-25,c("genotype", "pValue", "TT", "Symbol", "syn_annotation", "phenotype")]  %>% filter( pValue < 10^-4) %>% data.frame() %>%filter(syn_annotation == "Wobble_to_WCF") %>% head(10)

p2 <- ggplot(bin %>% filter(pValue < 10^-4) , aes(x = position, y = -log10(pValue + 1e-196), color = syn_annotation)) +
  geom_point(size = 0.5) + facet_grid( ~ syn_annotation) +
  pretty_plot(fontsize = 7) + scale_color_manual(values = jdb_palette("corona")) +
  theme(legend.position = "none") + labs(x = "Position in mtDNA genome", y = "-log10 p-value")
p2
cowplot::ggsave2(p2, file = "plots/association_binary_ryan2-updated.pdf", width = 4, height = 1.2)

cowplot::ggsave2(
  cowplot::plot_grid(p1, p2, ncol = 1, align = "v"),
  filename = "plots/stacked-association.pdf", 
  width = 4.3, height = 2.3)




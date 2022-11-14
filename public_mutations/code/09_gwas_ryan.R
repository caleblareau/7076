library(data.table)
library(dplyr)
library(BuenColors)


df <- fread('../data/gwas_ryan/mt_syn_binary.tsv') %>%
  filter(model == "Ryan2")

df %>% 
  filter(pValue < 1e-05) %>%
  group_by(genotype,syn_annotation) %>% summarize(n = n()) %>%
  mutate(Consequence = "synonymous_variant") %>% 
  count_me()

df %>% filter(pValue < 1e-5) %>% pull(genotype) %>% table() %>% sort()

df %>%
  filter(Codons == "ggG/ggA") %>%
  ggplot(aes(x = log10(oddsRatio), y = -log10(pValue))) +
  geom_point()

df %>% filter(Symbol == "MT-CO1")

dfq <- fread('../data/gwas_ryan/mt_syn_quant.tsv')
dfq %>%
  filter(Codons == "ggG/ggA")

ggplot(dfq, aes(x = syn_annotation, y = -log10(p)))  +
  geom_boxplot() +
  coord_cartesian(ylim = c(0,50))

dfq %>% 
  filter(pheno == "White blood cell (leukocyte) count") %>%
  ggplot(aes(x = syn_annotation, y = beta)) +
  geom_boxplot()

dfq %>% 
  filter(pheno == "Lactate") 

dfq %>%   
  filter(syn_annotation == "WCF_to_Wobble") %>%
  filter(p < 1e-08) %>% arrange((beta)) %>% data.frame() 
  

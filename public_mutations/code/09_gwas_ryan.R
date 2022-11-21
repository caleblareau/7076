library(data.table)
library(dplyr)
library(BuenColors)


df <- fread('../data/gwas_ryan/mt_syn_binary.tsv') %>%
  filter(model == "Ryan2")

df %>% filter(pValue < 10^-6) %>% pull(Symbol) %>% table()

df %>% filter(pValue < 10^-6) %>%
  group_by(Symbol, syn_annotation, genotype) %>% 
  summarize(count = n()) %>% ungroup() %>%
  group_by( syn_annotation) %>%
  summarize(count = n()) %>% data.frame()

df %>% 
  filter(pValue < 5e-05) %>%
  mutate(Consequence = "synonymous_variant") %>% 
  process_stats_syn_type()

df %>% filter(pValue < 1e-5) %>% 
  mutate(Consequence = "synonymous_variant") %>%
  count_me()

df %>%
  filter(Codons == "ggG/ggA") %>%
  ggplot(aes(x = log10(oddsRatio), y = -log10(pValue))) +
  geom_point()

df %>% filter(Symbol == "MT-CO1")

dfq <- fread('../data/gwas_ryan/mt_syn_quant.tsv')
dfq %>% filter(p < 10^-6) %>% 
  mutate(Consequence = "synonymous_variant") %>% 
  process_stats_syn_type()

dfq %>% filter(p < 10^-6) %>% group_by(Symbol, syn_annotation, geno) %>% 
  summarize(count = n()) %>% ungroup() %>%
  group_by(syn_annotation) %>%
  summarize(count = n()) %>% data.frame()

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


# Try kruskal wallis test
bdf <- data.frame(
  syn_annotation = c(df$syn_annotation, dfq$syn_annotation),
  p_value = c(df$pValue, dfq$p)
)
wilcox.test(bdf$p_value,bdf$syn_annotation=="WCF_to_WCF")

wilcox.test(bdf$p_value[bdf$syn_annotation=="WCF_to_WCF"],
            bdf$p_value[bdf$syn_annotation!="WCF_to_WCF"])


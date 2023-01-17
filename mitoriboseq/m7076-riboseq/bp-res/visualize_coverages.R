library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
process_tsv <- function(dt, xlength){
  avec <-  (dt %>%
              filter(V3 == "A") %>%
             # filter(V2 %in% c(33, 34)) %>% filter(V1 == 7046) %>% 
              mutate(chr = "chrM", start = V1, end = V1 + V2) %>%
              data.frame() %>%
              makeGRangesFromDataFrame() %>%
              coverage())[[1]] %>% as.numeric
  
  gvec <- (dt %>%
             filter(V3 == "G") %>%
             #filter(V2 %in% c(33, 34)) %>% filter(V1 == 7046) %>% 
             mutate(chr = "chrM", start = V1, end = V1 + V2 ) %>%
             data.frame() %>%
             makeGRangesFromDataFrame() %>%
             coverage())[[1]] %>% as.numeric
  
  data.frame(position = c(1:length(avec), 1:length(gvec)), coverage = c(avec, gvec), 
             letter = c(rep("A", length(avec)), rep("G", length(gvec)))) %>%
    filter(coverage > 0) %>%
    mutate(codon_color = case_when(
      position %in% c(7074,7075) ~ "zhetCodon",
      position %in% c(7076) ~ "het7076",
      position %in% (c(7074,7075,7076) +3) ~ "ACodon",
      TRUE ~ "other"
    ))
}

rnaseq <- rbind(fread("bams2/RNA_rep1_7076AG.tsv"), fread("bams2/RNA_rep2_7076AG.tsv"))
rbpseq <- rbind(fread("bams2/RBP_rep1_7076AG.tsv"), fread("bams2/RBP_rep2_7076AG.tsv"))

p1 <- process_tsv(rbpseq) %>%
  filter(position > 7050 & position < 7090) %>%
  ggplot(aes(x = position, y = coverage, fill = codon_color)) +
  geom_bar(color = "black", stat = "identity") +
  facet_wrap(~letter) +
  scale_fill_manual(values = c( "dodgerblue4", "firebrick", "lightgrey", "firebrick")) +
  scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 8) + theme(legend.position = "none")

cowplot::ggsave2(p1, file = "rbp_coverage.pdf", width = 3, height = 1.4)

14196/(9106+14196)
4474/(4474+5310)

process_tsv(rbpseq) %>% filter(position == 7076) 
process_tsv(rbpseq) %>% filter(position %in% c(7077:7079, 7080:7085)) %>% group_by(letter) %>%
  summarize(mean(coverage))


avec <-  (dt %>%
            filter(V3 == "A") %>%
           # filter(V2 %in% c(33, 34)) %>% # filter(V1 == 7046) %>% 
            mutate(chr = "chrM", start = V1 + V2 , end = V1 + V2) %>%
            data.frame() %>%
            makeGRangesFromDataFrame() %>%
            coverage())[[1]] %>% as.numeric

gvec <- (dt %>%
           filter(V3 == "G") %>%
           #filter(V2 %in% c(33, 34)) %>% # filter(V1 == 7046) %>% 
           mutate(chr = "chrM", start = V1 +V2, end = V1 +V2 ) %>%
           data.frame() %>%
           makeGRangesFromDataFrame() %>%
           coverage())[[1]] %>% as.numeric

data.frame(position = c(1:length(avec), 1:length(gvec)), coverage = c(avec, gvec), 
           letter = c(rep("A", length(avec)), rep("G", length(gvec)))) %>%
  filter(coverage > 0) %>%
  mutate(codon_color = case_when(
    position %in% c(7074,7075) ~ "zhetCodon",
    position %in% c(7076) ~ "het7076",
    position %in% (c(7074,7075,7076) +3) ~ "ACodon",
    TRUE ~ "other"
  )) %>% 
  mutate(x = 7076-position)

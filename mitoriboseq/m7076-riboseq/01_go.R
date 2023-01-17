library(data.table)
library(dplyr)
library(SummarizedExperiment)
library(Matrix)
library(ggbeeswarm)
library(BuenColors)

# Function that quickly computes the allele frequency matrix from a summarized experiment mgatk object
compute7076 <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  data.frame(
    nA = assays(SE)[[paste0("A_counts_rev")]][7076,] + assays(SE)[[paste0("A_counts_fw")]][7076,],
    nG = assays(SE)[[paste0("G_counts_rev")]][7076,] + assays(SE)[[paste0("G_counts_fw")]][7076,]
  ) %>%
    mutate(G_het = round((nG )/(nA + nG + 0.00001)*100, 1)) 
  
}

##### process
readRDS("round2-mitoribo.rds") %>% compute7076 %>%
  mutate(fraction = as.character(rep(c(0, 5:9),2))) -> ddf

fisher.test(
  matrix(c(sum(ddf$nA[c(1,7)]),sum(ddf$nG[c(1,7)]),
           sum(ddf$nA[c(2:6,8:12)]),sum(ddf$nG[c(2:6,8:12)])), nrow = 2)
) %>% str()


pX <- ddf %>%
  ggplot(aes(x = fraction, y = G_het)) +
  geom_bar(data = ddf %>% group_by(fraction) %>% summarize(mean=mean(G_het)), aes(y=mean), stat = "identity",
           color = "black", fill = "lightgrey", width = 0.7)+
  geom_point() + 
  coord_flip(ylim =c(40, 90)) + pretty_plot(fontsize = 7)  + L_border()
  
cowplot::ggsave2(pX, file = "barplot_newexp.pdf", width = 1.4, height = 1.5)



############################

# A / B are different donors; donor B is the relevant one here
readRDS("round1-mitoribo.rds") %>% compute7076

fisher.test(
  matrix(c(10212,6859, (301+208+45), (36+346+297)), nrow = 2)
) %>% str()

(36+346+297)/((301+208+45) + (36+346+297))
6859/(6859 + 10212)

mean_df <- data.frame(
  what = c("input", "riboprofiling"),
  value = c(40.18, 55.07)
)

(55.07-40.18)/40.18
rep_df <- data.frame(
  what = "riboprofiling",
  value = c(44.4, 53.8, 58.8)
)
pX <- ggplot(mean_df, aes(x = what, y = value)) +
  geom_bar(stat = "identity", fill = "lightgrey", color = "black", width = 0.7) +
  geom_quasirandom(data = rep_df, size = 0.5) +
  coord_flip(ylim =c(35, 65)) + pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Library", y = "m.7076A>G heteroplasmy (%)")

cowplot::ggsave2(pX, file = "histo_reps.pdf", 
                 width = 2.2, height = 1.6)



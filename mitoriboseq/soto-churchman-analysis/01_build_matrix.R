library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(pheatmap)
lf <- list.files(".", pattern = "*bedGraph*")

sapply(lf, function(x){
  import(x) %>%
  data.frame() -> df
  rep(df$score, df$width)
}) -> mat
colnames(mat) <- gsub(".bedGraph.gz", "", colnames(mat))
pheatmap(cor(mat))
mat[c(7075:7076),]

write.table(mat, file = "o.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

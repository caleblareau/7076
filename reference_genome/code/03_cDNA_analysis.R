library(data.table)
library(Biostrings)
library(stringr)
library(dplyr)

# Import and process mitochondria data
fa <- readDNAStringSet("../data/mt-sequences.fa")
cds_fa <- fa[width(fa) > 100]
cds_fa <- fa
mat <- str_split_fixed(str_split_fixed(str_split_fixed(names(cds_fa), "=", 4)[,4], "]", 2)[,1], "[.]", 2)
coords <- matrix(as.numeric(gsub("complement\\(", "", gsub(")", "", gsub("[.]", "", mat)))), ncol = 2)
lapply(1:dim(coords)[1], function(i){
  (coords[i,1]:coords[i,2])
}) %>% unlist() %>% unique() %>% length()

13854/16569
15358/16569

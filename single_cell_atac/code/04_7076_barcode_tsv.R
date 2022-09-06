library(data.table)
library(dplyr)
library(BuenColors)
"%ni%" <- Negate("%in%")

dates <- c(20201211, 20210115, 20210218, 20210329, 20210511)


files <- list.files("../data/mito_data/m7076/", full.names = TRUE)


the_expr <-  "../data/mito_data/m7076//ATS_PBMC_%s_C1_hg38_v20-7076.rds"

files <- lapply(dates, function(a) sprintf(the_expr,a)) %>% Reduce("c",.)
names(files) <- dates


## THIS FUNCTION ACCEPTS AN INDEX TO GIVE QUICK NAMES 
write_barcode_tsv <- function(i){
  
  
  ## Remove heteroplasmic cells (not 0 or 1) & filter for at least 10x coverage and assign allele value
  one <- readRDS(files[i]) %>% filter(af7076 %in% c(0,1), cov7076 > 10 ) %>% mutate(allele = if_else(af7076 ==0, "A","G"))
  
  
  ## Remove extra columns keeping only the two 
  one <- one %>% dplyr::select(cell, allele)
  
  new_out <- sprintf("../data/mito_data/mito_barcodes/ATS_PBMC_%s_allele_barcodes.tsv",names(files[i]))
  
  write.table(one, file = new_out, sep = "\t", quote = FALSE,col.names = FALSE,row.names = FALSE)
  
  
}


lapply(seq(length(files)),write_barcode_tsv)







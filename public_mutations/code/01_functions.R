library(dplyr)
library(BuenColors)
library(ggforce)
library(data.table)

annotate_synonomous <- function(df){
  sdf <- df %>% mutate(fp = Consequence == "synonymous_variant")
  variant_causes_Wobble_to_WCF <-  (grepl("\\*",sdf$Variant.tRNA) & grepl("near",sdf$Reference.tRNA)) & sdf$fp
  variant_causes_WCF_to_Wobble <-  (grepl("\\*",sdf$Reference.tRNA) & grepl("near",sdf$Variant.tRNA)) & sdf$fp
  variant_causes_Wobble_to_Wobble <-  (grepl("near",sdf$Variant.tRNA) & grepl("near",sdf$Reference.tRNA)) & sdf$fp
  variant_causes_WCF_to_WCF <- (grepl("\\*",sdf$Variant.tRNA) & grepl("\\*",sdf$Reference.tRNA)) & sdf$fp
  
  # Sanity check
  summary(as.numeric(variant_causes_Wobble_to_WCF) + 
            as.numeric(variant_causes_WCF_to_Wobble) + 
            as.numeric(variant_causes_Wobble_to_Wobble) + as.numeric(variant_causes_WCF_to_WCF))
  vec_annotate <- case_when(
    variant_causes_Wobble_to_WCF ~ "Wobble_to_WCF", 
    variant_causes_WCF_to_Wobble ~ "WCF_to_Wobble", 
    variant_causes_Wobble_to_Wobble ~ "Wobble_to_Wobble", 
    variant_causes_WCF_to_WCF ~ "WCF_to_WCF", 
    TRUE ~ "other"
  )
}

count_me <- function(df){
  df %>% 
    filter(Consequence == "synonymous_variant") %>%
    group_by(syn_annotation) %>%
    summarize(count = n()) %>%
    tidyr::complete(syn_annotation, fill = list(count = 0)) %>%
    mutate(perc = count/sum(count)*100)
}

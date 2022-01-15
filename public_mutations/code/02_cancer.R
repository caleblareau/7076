library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)

# Import Jacob's work
annotations <- fread("functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("m", annotations$Position, annotations$Reference, ">", annotations$Variant)

annotations[grepl("\\*",annotations$Reference.tRNA) & grepl("near",annotations$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>% data.frame() %>% dim()

# Import Phewas associations
assoc <- readxl::read_xlsx("MitoPhewas_associations.xlsx", 1) %>% data.frame()
assoc$mutation <- paste0("m", assoc$Position, toupper(assoc$Allele1),">",toupper(assoc$Allele2))

# Merge the two and look at FDR < 0.05 hits (how significance has been previously defined)
mdf <- merge(assoc, annotations, by = "mutation") %>% arrange(FDR)

kdf <- mdf[,c("Position.x","mutation","Symbol","Trait", "FDR", "Codons","Consequence", "Reference.tRNA", "Variant.tRNA", "Effect.size")] 
kdf05 <- kdf %>%
  filter(FDR < 0.05) 

kdf20 <- kdf %>%
  filter(FDR < 0.2) 


# Look at what recurrent phenotypes we have
kdf05[,c("mutation", "Consequence")]%>% group_by(mutation, Consequence) %>%
  summarize(count =n ()) %>% pull(Consequence) %>% table()

# Look at the patterns of mutations of interest
kdf20[grepl("\\*",kdf20$Reference.tRNA) & grepl("near",kdf20$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[grepl("\\*",kdf05$Reference.tRNA) & grepl("near",kdf05$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[grepl("\\*",kdf05$Variant.tRNA) & grepl("near",kdf05$Reference.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[grepl("near",kdf05$Reference.tRNA) & grepl("near",kdf05$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)


kdf05[grepl("\\*",kdf05$Reference.tRNA) & grepl("\\*",kdf05$Variant.tRNA)  ,] %>%
  filter(Consequence == "synonymous_variant") %>%
  arrange(Effect.size)

kdf05[,c("Trait", "Consequence")]%>% filter(Consequence == "synonymous_variant") %>%
  group_by(Trait) %>% summarize(count = n()) %>% arrange(desc(count))


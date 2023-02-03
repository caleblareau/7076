require(data.table)
require(dplyr)
require(stringr)
require(Biostrings)
rename <- dplyr::rename


## Read in MT tRNA anticodons
mt_trna_df <- readr::read_tsv("../../reference_genome/data/reference_tRNA_anticodons-wStop.tsv")
mt_trna_df <- mt_trna_df[mt_trna_df$AA.abr != "Stop",]
anti_codons_rna <- mt_trna_df$AA.anticodon %>% RNAStringSet()
match_trna_anticodons <- reverseComplement(anti_codons_rna) %>% DNAStringSet() 
mt_trna_df$AA.matchcodon <- as.character(match_trna_anticodons)


## JG 1/12/21: 
## The goal is to use the current data frame (all Codon) to match to variants. 
## I need to mark if the codon is a perfect match to a tRNA or the closest (Consider doing hamming distance, rn just mark if imperfect)
mt_trna_df$Codon.dna <- mt_trna_df$Codon %>% RNAStringSet() %>% DNAStringSet() %>% as.character()

mt_trna_df <- mt_trna_df %>% mutate(AA.matchbool = if_else(AA.matchcodon == Codon.dna,'*',"near"))



## Create a new index for only anticodons. 
#mt_trna_df$AA.idx <- sprintf("%s.%s",mt_trna_df$AA.abr,mt_trna_df$AA.anticodon)
mt_trna_df$AA.idx <- sprintf("%s.%s.%s",mt_trna_df$AA.matchbool,mt_trna_df$AA.abr,mt_trna_df$AA.anticodon)


## Now using Codon.dna to compare 
mt_trna_match_df <- mt_trna_df %>% select(Codon.dna,AA.idx)%>% 
  distinct() %>% dplyr::rename(Codon=Codon.dna, mt.tRNA = AA.idx) ## Get the 62 unique tRNA codon matches... THe other codons are stops I assume

## Now Add to ANnotations table 
functional_variant_df <- readr::read_tsv("mouse_mm10_synon_vars.tsv")
functional_variant_df <- functional_variant_df %>% tibble::rowid_to_column()


## SPlit it again by reference and variant. 
relevant_variant_df <- functional_variant_df %>% filter(Codons != "-")
## Process AA Data
aa_split <- relevant_variant_df$`Amino acids` %>% str_split_fixed("/",n=2)
empty_idx <- aa_split[,2] %>% grepl(pattern ="^$")


aa_split[,2][empty_idx] <- aa_split[,1][empty_idx] ## COpy the first (and only value) to make a complete matrix 

colnames(aa_split) <- c("Reference","Variant") %>% sprintf("%s.aa",.) 

## Now split codons and get the alt.position
codon_split <- relevant_variant_df$Codons %>% str_split_fixed("/",n=2)

#codon_split %>% melt(1:2) %>% pull(value) %>% unique()

#codon_split %>% as.data.frame() %>% transmute(across(everything(), ~grepl(.x,pattern =  "^A$")))

## Add the location of the AA change
codon_split <- cbind(codon_split,str_locate(string = codon_split[,1], pattern = "[A-Z]")[,"start"])

fake_mf_idx <- codon_split %>% as.data.frame() %>% transmute(across(everything(), ~grepl(.x,pattern =  "^A$"))) %>% Reduce("|",.) # %>% length()

#relevant_variant_df[fake_mf_idx,] %>% distinct(Consequence) %>% pull() %>% View()

## Format Data to uppercase and make it easy to match 
for (i in c(1,2) ){
  codon_split[,i] <- toupper(codon_split[,i])
  
}

col1_trna <- left_join(codon_split[,1,drop=F] %>% data.frame(Codon = .), mt_trna_match_df)
unmatched_codons <- col1_trna %>% filter(is.na(mt.tRNA))
col2_trna <- left_join(codon_split[,2,drop=F] %>% data.frame(Codon = .), mt_trna_match_df)


unmatched_codons <- rbind(unmatched_codons,col2_trna %>% filter(is.na(mt.tRNA)))
unmatched_codons <- unmatched_codons %>% distinct() %>% select(-mt.tRNA)

sprintf("Unmatched Codons: %s", paste(pull(unmatched_codons), collapse = ", ")) %>% noquote()

col1_trna <- col1_trna %>% rename(Reference.tRNA = mt.tRNA) %>% select(Reference.tRNA)
col2_trna <- col2_trna %>% rename(Variant.tRNA  = mt.tRNA) %>% select(Variant.tRNA)
relevant_variant_df <- data.frame(relevant_variant_df)
final_df <- cbind(relevant_variant_df[,c("X.Uploaded_variation","Location","Allele", "Codons", "Amino_acids")],col1_trna,col2_trna)
feat1 <- str_split_fixed(final_df$Reference.tRNA, "[.]", 3)[,1]
feat2 <- str_split_fixed(final_df$Variant.tRNA, "[.]", 3)[,1]

final_df$classify <- case_when(
  (feat1 == "near") & feat2 == "near" ~ "Wobble_to_Wobble",
  (feat1 == "*") & feat2 == "near" ~ "WCF_to_Wobble",
  (feat1 == "near") & feat2 == "*" ~ "Wobble_to_WCF",
  (feat1 == "*") & feat2 == "*" ~ "WCF_to_WCF"
)
final_df$Position <- as.numeric(str_split_fixed(final_df[["Location"]], "-", 2)[,2])
keep_vn <- c("Position", "X1", "X2", "Codons", "Amino_acids", "Reference.tRNA", "Variant.tRNA", "classify")
ffdf <- data.frame(
  final_df, str_split_fixed(str_split_fixed(final_df[["X.Uploaded_variation"]], "_", 3)[,3], "/", 2))[,c(keep_vn)]
colnames(ffdf) <- c("Position", "Ref", "Alt", "Codons", "Amino_acids", "Reference.tRNA", "Variant.tRNA", "classify")

table(ffdf$classify )/dim(ffdf)[1]

## Save new output for easy processing
write.table(ffdf, file = "mouse_variant_tRNA_anticodon_table.tsv", row.names = FALSE,
            col.names = TRUE, sep = "\t")

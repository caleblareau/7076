require(data.table)
require(dplyr)
require(stringr)
require(Biostrings)
rename <- dplyr::rename

## Read in Codon table and make clean
## Made into spooky function to record how the data was extracted. (copied the pdf into vim)
## NO LONGER USED!!! FOCUS ON THE NATRUE COMS PAPER WITH THE TRNA ANTICODONS. 
preprocess_jia_higgs <- function(){
  raw_codon_df <- read.table("Jia_Higgs_codontable.tsv",header = T)
  
  split_idx <- colnames(raw_codon_df) %>% str_split_fixed(pattern = "\\.",n=2)
  #split_idx[,2]
  #starting_idx <- split_idx[,2]  %>% grepl(pattern = "^$") ## Find the empty thign and replace it with 0
  split_idx[,2] <-  split_idx[,2]  %>% trimws() %>% gsub(pattern ="^$",replacement = "0") ## Put 0's in empty and split by this vector.
  
  
  #https://stackoverflow.com/questions/51959253/split-into-list-of-data-frames-by-column-index
  split_dfs <- split.default(raw_codon_df,f=split_idx[,2]) 
  new_colname_l <- split(split_idx[,1],split_idx[,2])
  
  
  split_dfs <- lapply(seq(length(split_dfs)) ,function(i){
    
    colnames(split_dfs[[i]]) <- new_colname_l[[i]]
    
    return(split_dfs[[i]])
    
  })
  
  long_df <- Reduce("rbind",split_dfs)
  
  write.table(long_df,"../data/annotations/Jia_Higgs_codontable_clean.tsv",row.names = F)
} ## End of silent preprocess function 


## Read in MT tRNA anticodons
mt_trna_df <- readr::read_tsv("nature_com_mt_trna_codons.txt")

keep_rows <- !is.na(mt_trna_df$AA.anticodon) ## Remove stop codons (THEY DONT HAVE ANTICODON MATCH)


mt_trna_df <- mt_trna_df[keep_rows,]


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

## This spits out 1's and 0's bc the hamming distances are so close and not helpful... 
# sapply(1:nrow(mt_trna_df) , function(i){
#   
#   q_v <- c(mt_trna_df[i,"AA.matchcodon"], mt_trna_df[i,"Codon.dna"]) %>% unlist()
#   d <- stringDist(q_v, method = "hamming")
#   return(d)
# }) -> hamming_v 




## JG 11/12/21: Obsolete :/ 
## Use this to append on AA.idx data to indicate the anticodon and AA.
#mt_trna_match_df <- mt_trna_df %>% select(AA.matchcodon,AA.idx)%>% 
  #distinct() %>% dplyr::rename(Codon=AA.matchcodon, mt.tRNA = AA.idx) ## Get the 22 unique tRNA codon matches... 


## Now using Codon.dna to compare 
mt_trna_match_df <- mt_trna_df %>% select(Codon.dna,AA.idx)%>% 
  distinct() %>% dplyr::rename(Codon=Codon.dna, mt.tRNA = AA.idx) ## Get the 62 unique tRNA codon matches... THe other codons are stops I assume



write.table(mt_trna_df, file = "../data/annotations/reference_tRNA_anticodons.tsv", row.names = F,sep = "\t")



## Now Add to ANnotations table 
functional_variant_df <- readr::read_tsv("../data/annotations/functional_variant_annotation_table.tsv")

functional_variant_df <- functional_variant_df %>% tibble::rowid_to_column()

## Look at entries wih non-'-' codons col and update to relevant data
#functional_variant_df %>% filter(Codons != "-") %>% head(10) %>% View() ## LETS DRIVE B)

## Looking at this I can split_fixed and if the second one is empty replace it with the first. 
#functional_variant_df %>% filter(Codons != "-") %>% count(`Amino acids`) %>% View()


## SPlit it again by reference and variant. 
relevant_variant_df <- functional_variant_df %>% filter(Codons != "-")
#rownames(relevant_variant_df) <- relevant_variant_df %>% mutate(rowname = sprintf("%i%s>%s",Position,Reference,Variant)) %>% pull(rowname)


#relevant_variant_df %>% filter(Position == 10760) %>% View()
## Make unique rownames to make back 
#rownames(relevant_variant_df) <- relevant_variant_df %>% mutate(rowname = sprintf("%i%s>%s",Position,Reference,Variant)) %>% pull(rowname)



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



## WTF is happening with this usage business 
## For the usage I will keep what it is currently? Just include the Reference.usage and Variant.usage
## I could try to include a normalized score? But what for idk... 
#col1_usage <- left_join(codon_split[,1,drop=F] %>% data.frame(Codon = .), long_df %>% select(Codon,Usage))

col1_trna <- left_join(codon_split[,1,drop=F] %>% data.frame(Codon = .), mt_trna_match_df)


unmatched_codons <- col1_trna %>% filter(is.na(mt.tRNA))


#col2_usage <- left_join(codon_split[,2,drop=F] %>% data.frame(Codon = .), long_df %>% select(Codon,Usage))
col2_trna <- left_join(codon_split[,2,drop=F] %>% data.frame(Codon = .), mt_trna_match_df)


unmatched_codons <- rbind(unmatched_codons,col2_trna %>% filter(is.na(mt.tRNA)))
unmatched_codons <- unmatched_codons %>% distinct() %>% select(-mt.tRNA)

sprintf("Unmatched Codons: %s", paste(pull(unmatched_codons), collapse = ", ")) %>% noquote()


#col1_usage <- col1_usage %>% rename(Reference.Usage = Usage) %>% select(Reference.Usage)
#col2_usage <- col2_usage %>% rename(Variant.Usage = Usage) %>% select(Variant.Usage)

col1_trna <- col1_trna %>% rename(Reference.tRNA = mt.tRNA) %>% select(Reference.tRNA)
col2_trna <- col2_trna %>% rename(Variant.tRNA  = mt.tRNA) %>% select(Variant.tRNA)

#relevant_variant_df <- cbind(relevant_variant_df,col1_usage,col2_usage)
relevant_variant_df <- cbind(relevant_variant_df,col1_trna,col2_trna)





## PARSE SIFT & POLY SCORE
#test <- relevant_variant_df$SIFT %>% stringr::str_extract_all("(?=\\().*?(?<=\\))",simplify = T)  %>% .[,1] %>% gsub(pattern = "\\(|\\)",replacement =  "")
extract_val <- function(x) x %>%  stringr::str_extract_all("(?=\\().*?(?<=\\))",simplify = T)  %>% .[,1] %>% gsub(pattern = "\\(|\\)",replacement =  "")
relevant_variant_df <- relevant_variant_df %>% mutate(SIFT.val= extract_val(SIFT), PolyPhen.val = extract_val(PolyPhen))


## Save new output for easy processing
write.table(relevant_variant_df, file = "../data/annotations/functional_variant_tRNA_anticodon_table.tsv", row.names = F,sep = "\t")


## Row_ids of weird codons
#weird_idxs <- codon_split[,1] %>% str_length() < 3

#relevant_variant_df$rowid[weird_idxs] %>% as.character() %>% dput() #%>% paste(collapse=", ")
#rownames(relevant_variant_df) <- relevant_variant_df$rowid

#relevant_variant_df[weird_idxs,]

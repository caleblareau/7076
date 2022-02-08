library(data.table)
library(Biostrings)
library(stringr)
library(dplyr)

# Import and process mitochondria data
fa <- readDNAStringSet("../data/mt-sequences.fa")
trna_fa <- fa[width(fa) < 100]

# Start parsing out the goodies
mat1 <- str_split_fixed(names(trna_fa), " ", 4)
trna_df <- data.frame(seq = unname(as.character(trna_fa)),
                      gene = gsub("\\]", "", gsub("\\[gene=","",mat1[,2]))
)

# Import reference data
reference_data <- fread("../data/coordinate_numbers.tsv") %>%
  mutate(gene = paste0("TRN", short1)) %>%
  mutate(width = Genome_Coord_2 - Genome_Coord_1)

reference_data$acpos1 = case_when(reference_data$Strand == "L" ~ as.numeric(reference_data$Genome_Coord_1 + reference_data$Pos_anticodon - 1),
                                                             TRUE ~ as.numeric(reference_data$Genome_Coord_1 + (reference_data$width - reference_data$Pos_anticodon) + 1))
reference_data$acpos2 = case_when(reference_data$Strand == "L" ~ reference_data$acpos1 + 1,
                                     TRUE ~ reference_data$acpos1 - 1)
reference_data$acpos3 = case_when(reference_data$Strand == "L" ~ reference_data$acpos1 + 2,
                                     TRUE ~ reference_data$acpos1 - 2)


# Verify anticodon from transcript
mdf <- merge(reference_data, trna_df, by = "gene")
mdf$parsed_anticodon <- substr(mdf$seq, mdf$Pos_anticodon, mdf$Pos_anticodon + 2)
mdf$parsed_anticodon == mdf$anticodon 

reference_data[,c("short3", "acpos1", "acpos2", "acpos3")]


# Verify positions from reference genome
ref_allele <- fread("../data/MT_refAllele.tsv")[["V2"]]
ref_allele_complement <- chartr("ATGC","TACG", ref_allele)

data.frame(
  first_pos_ref = ref_allele_complement[reference_data$acpos1],
  ac =reference_data$anticodon,
  last_pos_ref = ref_allele_complement[reference_data$acpos3],
  strand = reference_data$Strand
) %>% filter(strand == "H")

data.frame(
  first_pos_ref = ref_allele_complement[reference_data$acpos1],
  ac =reference_data$anticodon,
  last_pos_ref = ref_allele_complement[reference_data$acpos3],
  strand = reference_data$Strand
) %>% filter(strand == "L")

# Now export
reshape2::melt(reference_data[,c("short3", "acpos1", "acpos2", "acpos3")], id.vars = "short3") %>%
  arrange((value)) %>%
  write.table(quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", file = "../output/tRNA_anticodon_positions.tsv")


library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
library(VariantAnnotation)


# Now import tRNA annoations
trna_positions_wide <- read.table("../../reference_genome/data/mt-universal-trna-genome-coordinates.csv", sep = ",", header = TRUE) 
trna_positions_wide %>%
  reshape2::melt(id.vars = "Position") %>%
  filter(!is.na(value)) -> trna_positions

# Import cancer associations
impact <- fread("../../public_mutations/data/cancer/data_mutations_impact.txt.gz") %>%
  filter(Variant_Type == "SNP") %>%
  mutate(mutation = ShortVariantID)

tcga <- fread("../../public_mutations/data/cancer/data_mutations_tcga.txt.gz") %>%
  filter(Variant_Type == "SNP" & class == "somatic" & helix_freq < 0.1 & genbank_freq < 0.1) %>%
  mutate(mutation = ShortVariantID)

pcawg <- fread("../../public_mutations/data/cancer/data_mutations_pcawg.txt.gz") %>%
  filter(Variant_Type == "SNP") %>%
  mutate(mutation = paste0(Reference_Allele, Start_Position, Tumor_Seq_Allele2))

data.frame(
  mutation = c(impact$mutation, tcga$mutation, pcawg$mutation),
  position = c(impact$Start_Position, tcga$Start_Position, pcawg$Start_Position)
) %>%
  group_by(mutation,position) %>% summarize(count_tumors = n()) -> cancer_count_df

# Do analysis of tcga
trna_positions %>%
  filter(value %in% cancer_count_df$mutation)

mdf <- merge(cancer_count_df, trna_positions, by.x = "position", by.y = "value") %>%
  arrange((position))

ggplot(mdf, aes(x = as.numeric(as.character(Position)), y = variable, fill = count_tumors)) +
  geom_tile(color = "black")

mdf %>% group_by(Position) %>%
  summarize(sum = sum(count_tumors)) %>%
  ggplot(aes(x = as.numeric(as.character(Position)), y = sum)) +
  geom_point() + geom_line() +
  pretty_plot() + L_border() + labs(x = "Position in reference mt tRNA", y = "Total number of mutations observed")

mdf %>% group_by(Position) %>%
  summarize(sum = sum(count_tumors)) %>%
  arrange(as.numeric(as.character(Position))) %>% data.frame()

# Now look at mitomap
rna <- fread("../data/mitomap-RNA-mutations.txt")
positions_rna_all <- rna[complete.cases(rna),"V1"][[1]]
positions_rna_cfrm <- rna[rna$V8 == "Cfrm","V1"][[1]]

trna_positions$confirm <- trna_positions$value %in% positions_rna_cfrm
trna_positions$maybe <- trna_positions$value %in% positions_rna_all

trna_positions %>% group_by(Position) %>%
  summarize(cfrm = sum(confirm), all = sum(maybe)) %>%
  mutate(position = as.numeric(as.character(Position))) %>%
  reshape2::melt(id.vars = c("Position", "position")) %>% 
  ggplot(aes(x = position, y = value, color = variable)) +
  geom_point() + geom_line() + 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(x = "Position in reference tRNA", y = "# of mutation events", "mutation type") +
  pretty_plot() + L_border() + theme(legend.position = "bottom")

# Now look at gnomad
vcf <- readVcf( "../data/gnomad.genomes.v3.1.sites.chrM.vcf.bgz", verbose = FALSE )

simple_gnomad <- data.frame(rowRanges(vcf)@ranges,
                            n_homo = vcf@info$AC_hom) %>% filter(width == 1)
simple_gnomad$id <- paste0("m", simple_gnomad$start)
trna_positions$id <- paste0("m", trna_positions$value)

mdf <- merge(simple_gnomad, trna_positions, by = "id")
mdf2 <- mdf %>% group_by(Position, variable) %>%
  summarize(n_homo = sum(n_homo))
ggplot(mdf2, aes(x = as.numeric(as.character(Position)), y = variable, fill = log10(n_homo + 1), label = n_homo)) +
  geom_tile(color = "black") +geom_text() + scale_fill_gradientn(colors = jdb_palette("brewer_heat")) + pretty_plot() + L_border() +
  theme(legend.position = "bottom") + labs(x = "Position in reference tRNA", y = 'Specific tRNA', color = "log10 n_homoplasmic") 

mdf %>% group_by(Position) %>%
  summarize(sum = sum(n_homo)) %>%
  ggplot(aes(x = as.numeric(as.character(Position)), y = sum)) +
  geom_point() + geom_line() +
  pretty_plot() + L_border() + labs(x = "Position in reference mt tRNA", y = "Total number of mutations observed") +
  ggtitle("gnomad")

mdf %>% group_by(variable) %>%
  summarize(total = sum(n_homo)) %>% arrange(desc(total)) %>% data.frame()

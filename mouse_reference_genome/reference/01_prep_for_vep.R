library(data.table)
library(seqinr)
library(stringr)

# Go
ref = toupper(as.character(read.fasta("mm10.fasta")[["chrM"]]))
ref_var <- paste0(1:16299,"_", ref)
all_possible <- paste0(rep(1:16299, 4), "_", rep(c("A", "C", "G", "T"), each = 16299))

# do some munging
df <- data.frame(l = all_possible[!(all_possible %in% ref_var)])
mdf <- data.frame(df, str_split_fixed(df[[1]], "_", 2))
mdf$ref <- ref[as.numeric(mdf$X1)]

# Now format for bed
data.frame(
  chr = "MT",
  start = mdf$X1,
  end = mdf$X1,
  var = paste0(mdf$ref, "/", mdf$X2)
) %>% write.table(file = "all_vars.tsv",
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

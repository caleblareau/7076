library(data.table)
library(dplyr)
library(stringr)
library(BuenColors)

# Set up plotting preferences
order_df <- data.frame(
  rep(rep(c("T", "C", "A", "G"), each = 4),4),
  rep(c("T", "C", "A", "G"), each = 16),
  rep(c("T", "C", "A", "G"), 16))

collapse_df <- function(order_df){
  
  codon_order <- sapply(1:dim(order_df)[1], function(i){
    paste0(order_df[i,1], order_df[i,2], order_df[i,3])
  })
  codon_order
}
codon_order <- collapse_df(order_df)
codon_annotation_mito <- fread("../data/reference_tRNA_anticodons-wStop.tsv")
codon_annotation_nuc <- fread("../data/NUCLEAR-reference_tRNA_anticodons-wStop.txt")

# Compute for nuckear and mitochondria
raw <- fread("../data/o586358-Human_CDS.tsv")
table(raw$Organelle)
mito <- raw %>% filter(Organelle == "mitochondrion") %>% data.frame()
nuclear <- raw %>% filter(Organelle == "genomic")%>% data.frame()

# Compute fraction of 
mito_melt <- reshape2::melt(mito[,12:75])

wcf_codons <- codon_annotation_mito %>% filter(AA.matchbool == "*") %>%
  pull(Codon.dna)
mito_melt %>% pull(value) %>% sum() #3802 codons
mito_melt %>% filter(variable %in% wcf_codons) %>% pull(value) %>% sum() #1976 codons

(3802-1976)/3802 # 48% of codons require wobble effect

# Plotting stuff
compute_count_bias <- function(df, codon_annotation_df){
  
  # Merge with tRNA annotations
  mdf <- merge(
    data.frame(Number = colSums(data.matrix(df[,codon_order]))),
    codon_annotation_df, 
    by.x = "row.names", by.y = "Codon.dna") 
  
  mdf %>%
    group_by(AA.abr) %>%
    mutate(prop = Number/sum(Number), n_possible = n()) %>%
    ungroup() %>% mutate(bias = prop/(1/n_possible)) %>%
    mutate(log2_bias = log2(bias)) -> bdf 
  
  bdf$Codon <- factor(as.character(bdf$Row.names), codon_order)
  bdf <- bdf %>% arrange(Codon)
  bdf$first <- factor(substr(bdf$Codon,2,2), c("T", "C", "A", "G"))
  bdf$second <- rep(16:1, 4)
  bdf$bias_color <- case_when(
    bdf$log2_bias > 2 ~ 2, 
    bdf$log2_bias < -2 ~ -2,
    TRUE ~ bdf$log2_bias
  )
  bdf$label <- paste0(
    bdf$Codon, " ",  as.character(round(bdf$bias, 2))
  )
  bdf
}

nuc_bias <- compute_count_bias(nuclear, codon_annotation_nuc)
mito_bias <- compute_count_bias(mito, codon_annotation_mito)

### Do a viz of codon usage biases
aa_vec <- c("A","R","N","D","B","C","E","Q","Z","G","H","I","L","K","M","F","P","S","T","W","Y","V")
names(aa_vec) <- c("Ala","Arg","Asn","Asp","Asx","Cys","Glu","Gln","Glx","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val")

mito_bias %>%
  filter(AA.abr != "Stop") %>% 
  pull(log2_bias) %>%
  abs() %>% 
  mean()

nuc_bias %>%
  filter(AA.abr != "Stop") %>% 
  pull(log2_bias) %>%
  abs() %>% 
  mean()

mito_bias_plot <- mito_bias %>%
  filter(AA.abr != "Stop") %>%
  mutate(one = aa_vec[AA.abr]) %>%
  mutate(match = AA.matchbool != "near")

library(ggbeeswarm)
ggplot(mito_bias_plot, aes(x = match, y = bias, label = one)) +
  geom_boxplot() + geom_beeswarm(color = "dodgerblue") 

#----------
# Now visualize the results
xxtheme <-   theme(
  axis.line = element_blank(),
  axis.ticks.x = element_blank(),        ## <- this line
  axis.ticks.y = element_blank(),        ## <- this line
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(), 
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.title.y=element_blank(),
  axis.text.y=element_blank()
) 

x <- ifelse(mito_bias$AA.matchbool == "*", "*", "")
mito_bias$label <- paste0(mito_bias$label)

p1 <- ggplot(mito_bias, aes(x = first, y = second, label = label, fill = bias_color)) +
  geom_tile() + geom_text(size = 2) +
  scale_fill_gradient2(  low = ("firebrick"),
                         mid = "white",
                         high = ("dodgerblue4"),
                         midpoint = 0, limits = c(-2, 2)) + pretty_plot(fontsize = 6) + xxtheme  +
  theme(legend.position = "none")

p1

p2 <- ggplot(nuc_bias, aes(x = first, y = second, label = label, fill = bias_color)) +
  geom_tile() + geom_text(size = 2) +
  scale_fill_gradient2(  low = ("firebrick"),
                         mid = "white",
                         high = ("dodgerblue4"),
                         midpoint = 0, limits = c(-2, 2)) +
  pretty_plot(fontsize = 6) + xxtheme +
  theme(legend.position = "none") 

cowplot::ggsave2(cowplot::plot_grid(p1, p2), 
                 file = "../plots/codon_usage_bias_heatmap.pdf", 
                 width = 5, height = 3)


mito_bias %>% filter(AA.matchbool == "near") %>%
  arrange(desc(log2_bias))

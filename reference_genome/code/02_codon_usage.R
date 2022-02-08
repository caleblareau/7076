library(data.table)
library(dplyr)
library(stringr)

# Set up plotting preferences
order_df <- data.frame(rep(c("T", "C", "A", "G"), each = 16),
                       rep(rep(c("T", "C", "A", "G"), each = 4),4),
                       rep(c("T", "C", "A", "G"), 16))

collapse_df <- function(order_df){
  
  codon_order <- sapply(1:dim(order_df)[1], function(i){
    paste0(order_df[i,1], order_df[i,2], order_df[i,3])
  })
  codon_order
}
codon_order <- collapse_df(order_df)
codon_annotation <- fread("../data/reference_tRNA_anticodons-wStop.tsv")


# Compute for nuckear and mitochondria
raw <- fread("../data/o586358-Human_CDS.tsv")
table(raw$Organelle)
mito <- raw %>% filter(Organelle == "mitochondrion") %>% data.frame()
nuclear <- raw %>% filter(Organelle == "genomic")%>% data.frame()

# Plotting stuff
compute_count_bias <- function(df){
  
  # Merge with tRNA annotations
  mdf <- merge(
    data.frame(Number = colSums(data.matrix(df[,codon_order]))),
    codon_annotation, 
    by.x = "row.names", by.y = "Codon.dna") 
  
  mdf %>%
    group_by(AA.abr) %>%
    mutate(prop = Number/sum(Number), n_possible = n()) %>%
    ungroup() %>% mutate(bias = prop/(1/n_possible)) %>%
    mutate(log2_bias = log2(bias)) -> bdf 
  
  bdf$Codon <- factor(as.character(bdf$Row.names), codon_order)
  bdf <- bdf %>% arrange(Codon)
  bdf$first <- factor(substr(bdf$Codon,1,1), c("T", "C", "A", "G"))
  bdf$second <- rep(16:1, 4)
  bdf$bias_color <- case_when(
    bdf$log2_bias > 2 ~ 2, 
    bdf$log2_bias < -2 ~ -2,
    TRUE ~ bdf$log2_bias
  )
  bdf$label <- paste0(
    bdf$Codon, " ",  as.character(round(bdf$bias, 2)), " ", bdf$AA.abr
  )
  bdf
}

nuc_bias <- compute_count_bias(nuclear)
mito_bias <- compute_count_bias(mito)
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
mito_bias$label <- paste0(mito_bias$label, )x
p1 <- ggplot(mito_bias, aes(x = first, y = second, label = label, fill = bias_color)) +
  geom_tile() + geom_text() +
  scale_fill_gradient2(  low = ("firebrick"),
                         mid = "white",
                         high = ("dodgerblue4"),
                         midpoint = 0, limits = c(-2, 2)) + pretty_plot(fontsize = 6) + xxtheme  +
  theme(legend.position = "none") + ggtitle("Mitochondria")

p2 <- ggplot(nuc_bias, aes(x = first, y = second, label = label, fill = bias_color)) +
  geom_tile() + geom_text() +
  scale_fill_gradient2(  low = ("firebrick"),
                         mid = "white",
                         high = ("dodgerblue4"),
                         midpoint = 0, limits = c(-2, 2)) +
  pretty_plot(fontsize = 6) + xxtheme +
  theme(legend.position = "none") + ggtitle("Nuclear")

cowplot::plot_grid(p1, p2)


mito_bias %>% filter(AA.matchbool == "near") %>%
  arrange(desc(log2_bias))

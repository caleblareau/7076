library(data.table)
library(dplyr)
library(readxl)
library(BuenColors)
library(VariantAnnotation)
library(stringr)
source("01_functions.R")

hdb <- fread("../data/HelixMTdb_20200327.tsv") %>%
  filter(feature == "protein_coding_gene")

alleles <- str_split_fixed(hdb$alleles, ",", 2)
hdb$ref <- substr(alleles[,1],3,3)
hdb$alt <- substr(alleles[,2],2,2)

# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$id <- paste0(annotations$Reference,annotations$Position, annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

simple_helix <- data.frame(
  position = gsub("chrM:", "", hdb$locus),
  ref = hdb$ref, alt = hdb$alt,
  n_homo = hdb$counts_hom
) %>% mutate(id = paste0(ref, position, alt))

# Now merge
mdf <- merge(simple_helix, annotations, by = "id")
helix_homo <- mdf %>% filter(n_homo >=1)
helix_homo_common <- mdf %>% filter( n_homo >= 100)

### try something new
all <- mdf %>%
  filter(syn_annotation != "other") %>% 
  filter(syn_annotation == "Wobble_to_WCF") %>%
  group_by(Amino.acids) %>% 
  summarize(count = n())%>%
  mutate(prop = count / sum(count))

greater_100 <- mdf %>%
  filter(syn_annotation != "other") %>% 
  filter(syn_annotation == "Wobble_to_WCF") %>%
  filter(n_homo >= 100) %>% 
  group_by(Amino.acids) %>% 
  summarize(count = n()) %>%
  mutate(prophigh = count / sum(count))

annotations <- annotations %>%
  filter(syn_annotation != "other") %>% 
  filter(syn_annotation == "Wobble_to_WCF") %>%
  group_by(Amino.acids) %>% 
  summarize(count = n()) %>%
  mutate(prop_all_possible = count / sum(count))

merge_aa_df <- merge(merge(all, greater_100, by = "Amino.acids", all = TRUE),
                     annotations, by = "Amino.acids", all = TRUE)
merge_aa_df[is.na(merge_aa_df)] <- 0
merge_aa_df %>%
  ggplot(aes(x = prop, y = prophigh, label = Amino.acids)) + 
  geom_text()

library(viridis)

reshape2::melt(merge_aa_df[,c("Amino.acids","prop", "prop_all_possible")], id.vars = c("Amino.acids")) %>%
  ggplot(aes(x = variable, y = value*100, fill = Amino.acids)) + pretty_plot(fontsize = 8) +
  scale_y_continuous(expand = c(0,0)) + L_border() +
  geom_bar(stat = "identity", color = NA, width = 0.5)+ scale_fill_manual(values = rev(viridis(20))) +
  theme(legend.position = "none") + labs(x = "", y = "% vars") -> p1
cowplot::ggsave2(p1, file = "../output/stacked_bar_helix_aa.pdf", width = 1.2, height = 1.4)


chisq.test(merge_aa_df$count.y, merge_aa_df$count)

sum(merge_aa_df$count.y)
sum(merge_aa_df$count)

###########
# now do the standard enumeration
#######
process_stats_syn_type <- function(df){
  cm_a <- count_me(annotations)
  cm_t <- count_me(df)
  cm_t$all_perc <- cm_a$perc
  cm_t$pvalue <- c(prop.test(cm_t$count[1], sum(cm_t$count), p = cm_a$perc[1] / 100)$p.value,
                   prop.test(cm_t$count[2], sum(cm_t$count), p = cm_a$perc[2] / 100)$p.value,
                   prop.test(cm_t$count[3], sum(cm_t$count), p = cm_a$perc[3] / 100)$p.value,
                   prop.test(cm_t$count[4], sum(cm_t$count), p = cm_a$perc[4] / 100)$p.value)
  cm_t
  
}

process_stats_syn_type2 <- function(df){
  cm_t <- count_me(df)
  cm_t
  
}
process_stats_syn_type(helix_homo_common)
process_stats_syn_type(helix_homo)


# Make it continuous
count_me_helix_cutoff <- function(df, cutoff = 1){
  print(cutoff)
  df %>% 
    filter(n_homo >= cutoff) %>% 
    filter(Consequence == "synonymous_variant") %>%
    group_by(syn_annotation) %>%
    summarize(count = n()) %>%
    tidyr::complete(syn_annotation, fill = list(count = 0)) %>%
    mutate(perc = count/sum(count)*100) %>%
    mutate(null = process_stats_syn_type(annotations)[["perc"]]) %>%
    mutate(ratio = perc / null, diff = perc - null)
}

lapply(c(1:10, 15, 20, 30, 40, 50, 75, 100, 150, 200,  300, 400, 500), function(threshold){
  count_me_helix_cutoff(mdf, threshold) %>% mutate(n_carriers = threshold)
}) %>% rbindlist() -> all_df


Kendall::MannKendall(all_df %>% filter(syn_annotation == "Wobble_to_WCF") %>% pull(diff))

all_df %>% 
  ggplot(aes(x = n_carriers, y = (diff), color = syn_annotation)) + 
  geom_line() + 
  scale_x_log10() + 
  geom_point() +
  geom_hline(data = process_stats_syn_type(annotations) %>% mutate(diff = 0), aes(yintercept = diff), linetype = 2) +
  scale_y_continuous(limits = c(-20, 20)) +
  # facet_wrap(~syn_annotation, nrow = 1) + # scales = "free_y", 
  pretty_plot(fontsize = 7) +scale_color_manual(values = jdb_palette("corona"))+
  labs(x = "Minimum # of carriers of homoplasmic variant in Helixdb",
       y = "% of synonymous variants") + theme(legend.position = "none") -> p1
p1
cowplot::ggsave2(p1, file = "../output/trends_helix_cutoff_constantY_nofacet.pdf", width = 1.6, height = 1.5)


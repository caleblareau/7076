library(data.table)
library(dplyr)
library(BuenColors)
library(tidyverse)

####
# ´import nature genetics assoc
# Import Jacob's work
annotations <- fread("../data/functional_variant_tRNA_anticodon_table.tsv") %>% data.frame()
annotations$mutation <- paste0("M-", annotations$Position, "-",annotations$Reference, "-", annotations$Variant)
annotations$syn_annotation <- annotate_synonomous(annotations)

# Import NG paper associations
#W e next compared whether our test statistics correlated with those from a previously published array-based GWAS focused on MT-DNA variants
# in the UKB by Yonova-Doing et al. (2021) Note that this GWAS was performed on the full 500K cohort, compared to the 100K whole-genomes used here.
# We focused on replication of the quantitative traits since we’re better powered for those (Table S18). Due to differences in nomenclature,
# I manually annotated their supp table with phenotype labels that are consistent with ours (included under data/natgen_gwas/quant_assoc.tsv.
# There were three phenotypes they tested that we didn’t: eGFRcrcy, eGFRcr, and parent’s age.

# Helrpe
ReadQuantGWAS <- function(gwas.quant.fh, tested.vars.fh) {
  # Correct beta flips in the NG GWAS paper
  # Tested.variants: the variants we included in our ExWAS 
  gwas.quant <- fread(gwas.quant.fh)
  all.tested <- fread(tested.vars.fh) %>% 
    separate(variant, into = c("chr", "pos", "ref", "alt"), remove = F)
  
  # Get ref allele 
  ref.alleles <- all.tested %>% 
    mutate(pos = as.integer(pos))
  
  gwas.quant <- gwas.quant %>% 
    left_join(ref.alleles, by = c(Position = "pos"))
  
  # Flip beta if ref/alt are flipped
  gwas.quant <- gwas.quant %>% 
    mutate(alt = ifelse(`Effect Allele` != ref, `Effect Allele`, `Other allele`)) %>% 
    mutate(geno = paste0("M-", Position, "-", ref, "-", alt)) %>% 
    mutate(beta = ifelse(`Other allele` != ref, -1*beta, beta)) # Flip the beta
  
  return(gwas.quant)
}
gwas.quant <- ReadQuantGWAS("quant_assoc.txt", "all_tested_mtvars.csv")
wgs.quant <- fread("MT_quant_SYN_p0001.tsv")
gwas.quant <- gwas.quant %>% 
  left_join(wgs.quant, 
            by = c(geno = "geno", Pheno = "phenotype"), 
            suffix = c(".gwas", ".wgs")) %>% 
  filter(!is.na(Pheno))

gwas.sig <- gwas.quant %>% 
  filter(P < 5e-5)

library(ggpubr)
gwas.sig %>% 
  ggplot(aes(x=beta.gwas, y=beta.wgs)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("GWAS beta") + 
  ylab("WGS beta") +
  scale_x_continuous(limits = c(-0.05, 0.06)) + 
  scale_y_continuous(limits = c(-0.1, 0.1)) + 
  pretty_plot(fontsize = 7) + L_border() -> p1
cowplot::ggsave2(p1, file = "plots/effect_comparison.pdf", width = 1.8, height = 1.8)

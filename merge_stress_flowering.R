# set working directory
setwd("/Volumes/Franks_lab/johnson/results/enrichment/")

library(tidyverse)
library(magrittr)

#load stress genes and clean up names
stress_df <- read_csv("stress_genes.csv") %>% 
  set_colnames(
    c("idx", "tair_gene", "description", 
      "responsive_tf", "stress_elements", "chr")
    ) %>% 
  mutate(effect = "stress")

#load flowing time genes and clean up names
flowering_df <- read_csv("flowering_genes.csv") %>% 
  rename("tair_gene" = "Gene_details") %>% 
  mutate(effect = "flowering")


#load ortholog keys
ortho_df <- read_tsv("data/Athal_Brapa_genePairs.txt") %>% 
  set_colnames(c("tair_gene", "brad_gene"))


#MERGE!!!!
flowering_ortho <- inner_join(ortho_df, flowering_df, by = "tair_gene") %>% 
  select(tair_gene, brad_gene, effect)

stress_ortho <- inner_join(ortho_df, stress_df, by = "tair_gene") %>% 
  select(tair_gene, brad_gene, effect)

#stripped down df, with effect and gene pairs
effect_df <- full_join(flowering_ortho, stress_ortho) %>% distinct()

#get lifted exons with v1 gene names
gene_df <- 
  read_tsv(
    "data/Brapa_gene_v1tov3_lifted.bed", 
    col_names = c("chrom", "start", "end", "brad_gene")
  ) %>% 
  mutate(start = start - 1)

full_brad_df <- 
  inner_join(gene_df, effect_df, by="brad_gene") %>% 
  arrange(chrom, start, end)

filter(full_brad_df, effect == "stress")
filter(full_brad_df, effect == "flowering")

write_tsv(full_brad_df, "flowering_and_stress_genes.tsv", col_names = FALSE)

#14% overlap - seems ok - shoulder shrug. 
f_df <- split(full_brad_df, full_brad_df$effect)$flowering
s_df <- split(full_brad_df, full_brad_df$effect)$stress
mean(f_df$tair_gene %in% s_df$tair_gene)


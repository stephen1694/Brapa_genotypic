# set working directory
setwd("/Volumes/Franks_lab/johnson/results/enrichment/")

library(rvest)
library(tidyverse)
library(magrittr)

ft_page <- "http://www.phytosystems.ulg.ac.be/florid/databases/gene_list/flowering"
ft_html <- read_html(ft_page)

table_all <- ft_html %>% 
  rvest::html_nodes('body') %>% 
  xml2::xml_find_all(".//td | .//th") %>% 
  rvest::html_text()

table_ncol <- 9

ft_df <- 
1:table_ncol %>% 
  map_dfc(~{
    table_vec <- table_all[seq(.x, length(table_all) - table_ncol, by = table_ncol)]
    tibble(table_vec[-1]) %>% set_colnames(table_vec[1])
})

col_names_cleaned <- colnames(ft_df) %>% str_replace_all(" ", "_") %>% str_remove_all(pattern = "\t|\n")
ft_df <- set_colnames(ft_df, col_names_cleaned)

#downloaded from http://brassicadb.cn/#/FlowerGene/
brad_df <- readxl::read_excel("Flowering-related_genes_in_B.rapa.xlsx")
brad_names_cleaned <- colnames(brad_df) %>% str_replace_all(" ", "_") %>% str_remove_all("\\(|\\)")
brad_df <- set_colnames(brad_df, brad_names_cleaned)

brad_at_df <- 
full_join(ft_df, brad_df, by = c("Gene_details" = "At_ortholog"))
write_csv(brad_at_df, "flowering_genes.csv")




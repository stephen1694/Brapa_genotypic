# set working directory
setwd("/Volumes/Franks_lab/johnson/results/enrichment/")

library(rvest)
library(tidyverse)
library(magrittr)

full_at_df <- 
  1:5 %>% map_df(~{
    at_page <- str_glue("http://caps.ncbs.res.in/cgi-bin/mini/databases/stifdb2/fetch_gene_list_by_chromosome.pl?orgid=aratha&chr=AT{.x}")
    at_html <- read_html(at_page)
    
    at_table <- at_html %>% 
      rvest::html_nodes('.geneDetails_table') %>%
      xml2::xml_find_all(".//td | .//th") %>% 
      rvest::html_text()
    
    table_ncol <- 5
    
    at_df <- 
      1:table_ncol %>% 
      map_dfc(~{
        table_vec <- at_table[seq(.x, length(at_table) - table_ncol, by = table_ncol)]
        tibble(table_vec[-1]) %>% set_colnames(table_vec[1])
      })
    at_df$`TAIR ID` <- str_remove_all(at_df$`TAIR ID`, "[\n ]")
    at_df$chr <- .x
    at_df
  })

write_csv(full_at_df, "stress_genes.csv")

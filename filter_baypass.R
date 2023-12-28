# set working directory
setwd("/Volumes/Franks_lab/johnson/results/baypass/")

library(dplyr)
library(data.table)

filter_baypass <- function(file_in, vep_out, baypass_out, thin_n){
  #grab function argument
  
  #read files in
  vep <- fread(file_in, header = F, 
               showProgress = T)
  
  colnames(vep)[1:5] <- c("chrom", "pos", "pos2", "ref_alt", "freq")
  
  #bayp <- fread(input = baypass_in, header = F, showProgress = T)

  #merge columns
    #add:
  #column with chromosome name and position
  #column of distance between adjacent rows
  #column with row index
  vep_bayp <- as_tibble(vep) %>% 
    mutate(chrompos = paste0(chrom, pos),
           pos_dist = c(pos[2:length(pos)], NA) - c(pos[1:length(pos)]),
           idx = 1:n()
    )
  
  #vep_bayp$pos_dist
  
  #create vector that will hold repeated chrompos
  idx_rep <- vep_bayp$chrompos
  
  #if the distance is smaller than cutoff, keep previous chrompos as current
  for(i in 1:(length(idx_rep)-1)){
    if(vep_bayp$pos_dist[i] < thin_n){
      idx_rep[i+1] <- idx_rep[i]
    }
  }
  
  #overwrite vep_bayp
    #if chrompos is shared 
  vep_bayp %<>% 
    mutate(idx_repeated = idx_rep) %>% 
    mutate(freq_var = (1-freq)*freq) %>% 
    group_by(idx_repeated) %>% 
    filter(freq_var == max(freq_var)) %>% 
    ungroup()
  
  
  bay_filt <- vep_bayp %>% select(starts_with("V"))
  vep_filt <- vep_bayp %>% select(c("chrom", "pos", "pos2", "ref_alt", "freq"))
  
  fwrite(bay_filt, 
         file = baypass_out, 
         quote = F, 
         col.names = F, 
         sep = " "
  )
  
  fwrite(vep_filt, 
         file = vep_out, 
         quote = F, 
         col.names = F, 
         sep = " "
  )
  
}

filter_baypass(
  thin_n = 100,
  file_in = "anc_vs_desc_vep_baypass_r0.95_d50_L15_M200_q30_Q30.txt",
  vep_out = "anc_vs_desc_vep_r0.95_d50_L15_M200_q30_Q30_thin100.txt", 
  baypass_out ="anc_vs_desc_baypass_r0.95_d50_L15_M200_q30_Q30_thin100.txt"
)


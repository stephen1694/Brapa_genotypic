### Effective population size estimate  from differences in reproduction (Table S6)

# set working directory
setwd("/Volumes/Franks_lab/johnson/results")

# load packages
library(dplyr) # used to group rows by column values (group_by), filter values (filter), count values (tally), calculate statistics (summarize), and for pipe (%>%)

# set seed
set.seed(19150)

# loading in data file with traits including fitness as data frame
test_gen_df <-read.csv("Johnson_et_al_2022_Evolution_individuals_data.csv", header = T, stringsAsFactors = FALSE)	
str(test_gen_df)

# subset of 2 ancestor pops of 35 under drought treat to explore
drought <- test_gen_df[1:280,c(9,31)]
drought[is.na(drought)] <- 0

# subset of 2 ancestor pops of 35 under watered treat to explore
watered <- test_gen_df[281:560,c(9,31)]
watered[is.na(watered)] <- 0

# create data frame (Ne to store final Ne of bootstraps, fitness to store info within each iteration)
Ne <- data.frame(matrix(ncol = 1, nrow = 100))
colnames(Ne) <- "N_eff"
fitness <- data.frame(matrix(ncol = 1, nrow = 73))
colnames(fitness) <- "Fitness"

# function to estimate Ne from reproductive variance
new.function <- function(group){
  for (k in 1:1000){
    # randomly sample 73 individuals from data frame without replacement
    fit_p <- as.data.frame(table(sample(nrow(group), size=73, replace = FALSE)))
    var <- as.numeric(fit_p$Var1)
    
    # add the relative fitness of each individual selected (list in var) to the second column of fit_p
    for (i in var){
      fit_p[i,2] <- group[i,2]
    }
    
    # replacing na (individuals selected that did not reproduce) with zeros and name columns
    fit_p[is.na(fit_p)] <- 0
    colnames(fit_p) <- c("row", "Rel_seed_mass")
    
    # calculate total fitness of 73 sampled individuals and calculate relative fitness
    total <- sum(fit_p$Rel_seed_mass)
    for (z in 1:73){
      fit_p[1,2] <- fit_p[z,2] * (73 / total)
    }
    
    # sample the list of selected individuals weighted bythe new relative fitness calculated, with  replacement
    fit <- as.data.frame(table(sample(seq_len(nrow(fit_p)), size=73, prob = fit_p$Rel_seed_mass, replace = TRUE)))
    colnames(fit) <- c("Individual", "Fitness")
    fit$Individual <- as.character(fit$Individual)
    fit$Individual <- as.numeric(fit$Individual)
    
    # add offspring number to fitness data frame, replace NA with zeros
    for (i in 1:nrow(fit)){
      j <- fit[i,1]
      fitness[j,1] <- fit[i,2]}
    fitness[is.na(fitness)] <- 0
    
    # calculate variance in relative offspring produced
    var_df <- as.data.frame(fitness %>% dplyr::summarize(Mean = mean(Fitness), sd = sd(Fitness)))
    var_df[,2] <- var_df[,2] ^2
    colnames(var_df) <- c("Mean", "Variance")
    
    # calculate estimate of Ne
    Ne[k,1] <- (4*73) / (2+var_df[1,2])
  }
  
  # sort and report Ne
  Ne <- Ne %>% arrange(N_eff)
  my_list <- list("mean" = mean(Ne[,1]), "2.5%" = Ne[26,1], "97.5%" = Ne[975,1], "range" = range(Ne))
  return(my_list) 
}

# calculate for drought treatment (reported in Table S6)
new.function(group=drought)
# mean of 83.0 with 95% CI of 70.5-94.4 and range of 64.0-100.1

# calculate for watered treatment (reported in Table S6)
new.function(group=watered)
# mean of 106.8 with 95% CI of 93.9-119.3 and range of 83.4-129.4

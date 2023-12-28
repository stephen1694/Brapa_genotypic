### Tajimas pi nucleotide diversity (Figure S1 and Table S5)

# set working directory
setwd("/Volumes/Franks_lab/johnson/results/pi")

# loading libraries
library(vroom)
library(tidyverse)
library(reshape2)

# making list of chromosomes to filter scaffolds
chr_list <- c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10")

## part 4A: Tajima's pi for by regime (Figure S1)

# reading in ancestor pi data to ancestor_pi, filtering out scaffolds
ancestor_pi <- as.data.frame(vroom("Ancestor.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "A_pi")))
ancestor_pi <- ancestor_pi %>% filter(chromosome %in% chr_list) %>% mutate(A_pi = as.numeric(A_pi))

# preparing ancestor_pi: group by chromosome, computes chromosome size, calculates cululative start position of each chromosome, adds this info (in tot column) to dataset, calculates cumulative position of each SNP and adds to position_cum column
ancestor_pi <- ancestor_pi %>% 
  group_by(chromosome) %>% 
  summarise(chr_len=max(position)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(ancestor_pi, ., by=c("chromosome"="chromosome")) %>%
  mutate(position_cum=position+tot)
str(ancestor_pi)

# reading in pi data to drought_pi, filtering out scaffolds
drought_pi <- as.data.frame(vroom("Drought.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "DD_pi")))
drought_pi <- drought_pi %>% filter(chromosome %in% chr_list) %>% mutate(DD_pi = as.numeric(DD_pi))

# preparing drought_pi: group by chromosome, computes chromosome size, calculates cululative start position of each chromosome, adds this info (in tot column) to dataset, calculates cumulative position of each SNP and adds to position_cum column
drought_pi <- drought_pi %>% 
  group_by(chromosome) %>% 
  summarise(chr_len=max(position)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(drought_pi, ., by=c("chromosome"="chromosome")) %>%
  mutate(position_cum=position+tot)
str(drought_pi)

# reading in pi data to watered_pi, filtering out scaffolds
watered_pi <- as.data.frame(vroom("Watered.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "WD_pi")))
watered_pi <- watered_pi %>% filter(chromosome %in% chr_list) %>% mutate(WD_pi = as.numeric(WD_pi))

# preparing watered_pi: group by chromosome, computes chromosome size, calculates cululative start position of each chromosome, adds this info (in tot column) to dataset, calculates cumulative position of each SNP and adds to position_cum column
watered_pi <- watered_pi %>% 
  group_by(chromosome) %>% 
  summarise(chr_len=max(position)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(watered_pi, ., by=c("chromosome"="chromosome")) %>%
  mutate(position_cum=position+tot)
str(watered_pi)

# create one df with all pi values
merged_pi <- cbind(ancestor_pi[1],ancestor_pi[7],ancestor_pi[5],drought_pi[5],watered_pi[5])
str(merged_pi)

# checking means
mean(merged_pi$A_pi, na.rm=TRUE)
#0.01051094 for ancestor regime
mean(merged_pi$DD_pi, na.rm=TRUE)
#0.0103569 for drought regime
mean(merged_pi$WD_pi, na.rm=TRUE)
#0.01039544 for watered regime

# checking sd
sd(merged_pi$A_pi, na.rm=TRUE)
#0.003962348
sd(merged_pi$DD_pi, na.rm=TRUE)
#0.004036341
sd(merged_pi$WD_pi, na.rm=TRUE)
#0.003906866

# make new df with cols 2-5 and name
merged_pi2 <- merged_pi[,2:5]
colnames(merged_pi2) <- c("position", "Ancestor", "Drought", "Watered")

# melting so that all pi measurements are in 1 column named pi and regime information is provided in variable column
melted_pi_df <- melt(merged_pi2, id.vars = 'position')
colnames(melted_pi_df) <- c("position", "variable", "pi")

# printing Figure S1
FigS1 <- ggplot(melted_pi_df, aes(x=position, y=pi, fill=variable)) + geom_line() + facet_grid(variable ~ .) + scale_x_continuous(label = axisdf$chromosome, breaks= axisdf$center, expand = c(0,0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.045)) + theme_classic() + theme(legend.position = "n")
FigS1
ggplot2::ggsave(filename = "FigS1.pdf", plot = FigS1, width = 18.4, height = 18.4, units = "cm", device='pdf', dpi=300)


## part 4B: Tajimas pi for each separate drought pop (Table S5)

# reading in drought rep pi data to drought_pi, filtering out scaffolds
drought1_pi <- as.data.frame(vroom("Library1.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D1_pi")))
drought1_pi <- drought1_pi %>% filter(chromosome %in% chr_list) %>% mutate(D1_pi = as.numeric(D1_pi))
drought2_pi <- as.data.frame(vroom("Library2.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D2_pi")))
drought2_pi <- drought2_pi %>% filter(chromosome %in% chr_list) %>% mutate(D2_pi = as.numeric(D2_pi))
drought3_pi <- as.data.frame(vroom("Library3.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D3_pi")))
drought3_pi <- drought3_pi %>% filter(chromosome %in% chr_list) %>% mutate(D3_pi = as.numeric(D3_pi))
drought4_pi <- as.data.frame(vroom("Library4.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D4_pi")))
drought4_pi <- drought4_pi %>% filter(chromosome %in% chr_list) %>% mutate(D4_pi = as.numeric(D4_pi))
drought5_pi <- as.data.frame(vroom("Library5.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D5_pi")))
drought5_pi <- drought5_pi %>% filter(chromosome %in% chr_list) %>% mutate(D5_pi = as.numeric(D5_pi))
drought6_pi <- as.data.frame(vroom("Library6.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D6_pi")))
drought6_pi <- drought6_pi %>% filter(chromosome %in% chr_list) %>% mutate(D6_pi = as.numeric(D6_pi))
drought7_pi <- as.data.frame(vroom("Library7.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D7_pi")))
drought7_pi <- drought7_pi %>% filter(chromosome %in% chr_list) %>% mutate(D7_pi = as.numeric(D7_pi))
drought8_pi <- as.data.frame(vroom("Library8.pi", col_names = c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "D8_pi")))
drought8_pi <- drought8_pi %>% filter(chromosome %in% chr_list) %>% mutate(D8_pi = as.numeric(D8_pi))

# create empty df for function return
pi <- data.frame(matrix(NA, nrow = 1, ncol = 2))
colnames(pi) <- c("mean", "sd")

# function to calculate mean and sd of pi
pi_calc <- function(drought_pop){
  # preparing drought1_pi: group by chromosome, computes chromosome size, calculates cululative start position of each chromosome, adds this info (in tot column) to dataset, calculates cumulative position of each SNP and adds to position_cum column
  drought_pop_pi <- drought_pop %>% 
    group_by(chromosome) %>% 
    summarise(chr_len=max(position)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(drought_pop, ., by=c("chromosome"="chromosome")) %>%
    mutate(position_cum=position+tot)
  pi[1,1] <- mean(drought_pop_pi[,5], na.rm=TRUE)
  pi[1,2] <- sd(drought_pop_pi[,5], na.rm=TRUE)
  return(pi)
}

# call function for each drought library (reported in Table S5)
pi_calc(drought1_pi)
pi_calc(drought2_pi)
pi_calc(drought3_pi)
pi_calc(drought4_pi)
pi_calc(drought5_pi)
pi_calc(drought6_pi)
pi_calc(drought7_pi)
pi_calc(drought8_pi)

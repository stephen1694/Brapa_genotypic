### Manhattan plots (Figure 2)

# load required packages
library(tidyverse)
library(vroom)
library(magrittr)
library(gridExtra)

# making list of chromosomes
chr_list <- c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10")


## part 2A: manhattan plot for bayesiam model results ran with BayPass

# set working directory
setwd("/Volumes/Franks_lab/johnson/results/baypass/")

# read in data files 
betai_df <- read.table("anc_vs_desc_aux_summary_betai.out", header = T)
sites_df <- read.table("anc_vs_desc_vep_r0.95_d50_L15_M200_q30_Q30_thin100.txt") %>% set_colnames(c("chrom", "pos", "pos2", "ref_alt", "freq")) %>% mutate(MRK = 1:n())

# join and inspect 
full_df <- full_join(sites_df, betai_df, by = "MRK")
str(full_df)

# saved this table
write.table(full_df, file = "anc_vs_desc_full_df.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

# setting up df
filt_df <- full_df %>% 
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(full_df, ., by=c("chrom"="chrom")) %>%
  mutate(position_cum=pos+tot) %>% filter(BF.dB. >= 0) %>% filter(chrom %in% chr_list)
str(filt_df)

# setting x axis labels
axisdf <- as.data.frame(filt_df %>% group_by(chrom) %>% summarize(center=( max(position_cum) + min(position_cum) ) / 2 )) %>% filter(chrom %in% chr_list)

# count SNPs reaching sig threshholds from each covariate
filter(filt_df, BF.dB. >= 10 & COVARIABLE == 1) %>% count() #466
filter(filt_df, BF.dB. >= 20 & COVARIABLE == 1) %>% count() #51
filter(filt_df, BF.dB. >= 10 & COVARIABLE ==2) %>% count() #0
# none from cov 2 (drought vs watered), so will just plot cov 1 (ancestor vs descendants)

# preparing plot
cov1_plot <- ggplot(data = filter(filt_df, COVARIABLE == 1), aes(x=position_cum, y=BF.dB., color=as.factor(chrom))) +
  geom_point(alpha = 0.8, size = 1) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand =c(0,0)) +
  scale_color_manual(values = rep(c("black", "#276FBF"), unique(length(axisdf$chrom)))) +
  geom_hline(yintercept = 10, lty = 2) +
  geom_hline(yintercept = 20, lty = 2) +
  labs(x="Position", y="Bayes factor (dB)", tag ="A") +
  scale_y_continuous(expand = c(0, 0), limits = c(-14,54), breaks = c(-10,0,10,20,30,40,50)) + 
  theme_classic(base_size=12) +
  theme(legend.position = "n") 
# can inspect with print(cov1_plot)


## part 2B: manhattan plot for fisher exact tests (FETs) of drought libraries and ancestor regime

# set working directory
setwd("/Volumes/Franks_lab/johnson/results/FET/ancestor_regime_vs_drought_pops")

# read in data to brassica_df, remove characters from p value columns, convert the resulting p values in these columns to numeric, filter by chromosomes listed in chr_list to get rid of SNPs on scaffolds, and finally convert to data frame
brassica_df <- as.data.frame(vroom("Ancestor_drought_pops.fet"))
brassica_df <- brassica_df[,1:13]
colnames(brassica_df) <- c("chromosome", "position", "num_SNPs", "prop_SNPs_coverage", "min_coverage", "D1_vs_AH", "D2_vs_AH", "D3_vs_AH", "D4_vs_AH", "D5_vs_AH", "D6_vs_AH", "D7_vs_AH", "D8_vs_AH")

# removing unwanted characters and converting numbers to numeric
brassica_df <- brassica_df %>% mutate(D1_vs_AH = as.numeric(gsub("1:2=","",as.character(D1_vs_AH)))) %>% mutate(D2_vs_AH = as.numeric(gsub("1:3=","",as.character(D2_vs_AH)))) %>% mutate(D3_vs_AH = as.numeric(gsub("1:4=","",as.character(D3_vs_AH)))) %>% mutate(D4_vs_AH = as.numeric(gsub("1:5=","",as.character(D4_vs_AH)))) %>% mutate(D5_vs_AH = as.numeric(gsub("1:6=","",as.character(D5_vs_AH)))) %>% mutate(D6_vs_AH = as.numeric(gsub("1:7=","",as.character(D6_vs_AH)))) %>% mutate(D7_vs_AH = as.numeric(gsub("1:8=","",as.character(D7_vs_AH)))) %>% mutate(D8_vs_AH = as.numeric(gsub("1:9=","",as.character(D8_vs_AH)))) %>% filter(chromosome %in% chr_list)
str(brassica_df)

# preparing brassica_df: group by chromosome, compute chromosome size, calculate cumulative start position of each chromosome, add this info (in tot column) to dataset, calculate cumulative position of each SNP and add this with position_cum column
brassica_df <- brassica_df %>% 
  group_by(chromosome) %>% 
  summarise(chr_len=max(position)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(brassica_df, ., by=c("chromosome"="chromosome")) %>%
  mutate(position_cum=position+tot)

# calculating significance threshold based on # SNPs (bonferroni correction)
sig_cut <- -log10(0.05 / nrow(brassica_df))
sig_cut <- -log10(0.05 / (nrow(brassica_df) * 8))

# setting x axis labels
axisdf <- as.data.frame(brassica_df %>% group_by(chromosome) %>% summarize(center=( max(position_cum) + min(position_cum) ) / 2 ))

# settin ylim for graph
for (i in (6:13)){
  print(max(brassica_df[i]))}
ylim <- 35

# counting significant SNPs
for (i in (6:13)){
  column <- brassica_df[i]
  num <- sum(column >= sig_cut)
  prop <- (sum(column >= sig_cut) / nrow(brassica_df)) * 100
  mystring <- "SNPs, or a percentage of"
  print(paste(num, mystring, round(prop,2)))}
# "6435 SNPs, or a percentage of 0.14"
# "105 SNPs, or a percentage of 0"
# "310 SNPs, or a percentage of 0.01"
# "7126 SNPs, or a percentage of 0.16"
# "701 SNPs, or a percentage of 0.02"
# "261 SNPs, or a percentage of 0.01"
# "2998 SNPs, or a percentage of 0.07"
# "2967 SNPs, or a percentage of 0.07"

# combine all rows from columns 6-13 into 1 column for qq plot, manhattan plot analysis
new_df <- as.data.frame(rep(brassica_df$chromosome, times = 8))
new_df[,2] <- rep(brassica_df$position, times = 8)
new_df[,3] <- rep(brassica_df$position, times = 8)
new_df[,2] <- new_df[,2] -1
new_df[,4] <- c(rep("D1_vs_AH", times = 4479319), rep("D2_vs_AH", times = 4479319),rep("D3_vs_AH", times = 4479319), rep("D4_vs_AH", times = 4479319),rep("D5_vs_AH", times = 4479319), rep("D6_vs_AH", times = 4479319),rep("D7_vs_AH", times = 4479319), rep("D8_vs_AH", times = 4479319))
new_df[,5] <- c(brassica_df$D1_vs_AH, brassica_df$D2_vs_AH, brassica_df$D3_vs_AH, brassica_df$D4_vs_AH, brassica_df$D5_vs_AH, brassica_df$D6_vs_AH, brassica_df$D7_vs_AH, brassica_df$D8_vs_AH)
new_df[,6] <- rep(brassica_df$position_cum, times = 8)
colnames(new_df) <- c("chrom", "start", "end", "pop", "pvalue", "position_cum")
str(new_df)

# defining uniform_qq function
uniform_qq <- function(p_obs, err, log_p = TRUE){
  if(log_p){
    p_obs <- 10^(-p_obs)
  }
  qs <- (1:length(p_obs)-0.5) / length(p_obs)
  my_q <- qunif(qs)
  qq_df <- data.frame(
    theor =-log10(my_q), obs = -log10(quantile(p_obs, probs = qs))) %>%
    mutate(sq_diff = (theor - obs)^2) %>%
    filter(sq_diff > err)
  return(qq_df)
}

# calculating genomic inflation factor (lambda) for P values
chisq <- qchisq(10^-(new_df$pvalue), df=1, lower.tail = FALSE)
lambda <- median(chisq)/qchisq(0.5,1)
lambda
#1.62

# correcting chi square test values by lambda and calculating new -log10(P)
newchisq <- chisq/lambda
new_df$pvalue_new <- -log10(pchisq(newchisq, df=1, lower.tail = FALSE))
str(new_df)

# calculating genomic inflation factor (lambda) on corrected P vals
chisq1 <- qchisq(10^-(new_df$pvalue_new), df=1, lower.tail = FALSE)
lambda1 <- median(chisq1)/qchisq(0.5,1)
lambda1
#1.0

# preparing manhattan plot of GC drought descendant pops vs ancestor Regime SNP shifts 
DDvsAR_plot <- ggplot(data = filter(new_df, pvalue_new > 5), aes(x=position_cum, y=pvalue_new, color=as.factor(chrom))) +
  geom_point(alpha = 0.8, size=1) +
  scale_x_continuous(label = axisdf$chromosome, breaks= axisdf$center, expand = c(0,0)) +
  scale_color_manual(values = rep(c("black", "#276FBF"), unique(length(axisdf$chromosome)))) +
  geom_hline(yintercept = sig_cut, lty = 2) +
  labs(x="Position", y="-log10 P-value", tag ="B") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  theme_classic(base_size=12) +
  theme(legend.position = "n")
# can inspect with print(DDvsAR_plot)

## list of sig SNPs
SNP_diff_GC <- new_df %>% filter(new_df$pvalue_new >= sig_cut)
str(SNP_diff_GC) #913 differentiated SNPs
length(unique(SNP_diff_GC[["position_cum"]])) #911 unique sites

# save df from R to .txt file in FET directory
write.table(new_df, file = "anc_vs_drought_full_df.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)


## part 2C: printing Figure 2 (Baypass manhattan plot in panel A and fishers exact test manhattan plot in panel B)
Fig2 <- gridExtra::grid.arrange(cov1_plot, DDvsAR_plot, nrow = 2)
Fig2
ggplot2::ggsave(filename = "Fig2.pdf", plot = Fig2, width = 18.4, height = 18.4, units = "cm", device='pdf', dpi=300)

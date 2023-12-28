### SNP frequency PCA (Figure 1)

# load required packages
library(vroom)
library(tidyverse)
library(ggfortify)

# set working directory
setwd("/Volumes/Franks_lab/johnson/results/snp_diff/")

# load in allele counts
brassica_SNPs <- as.data.frame(vroom("Population_level_rc", col_names = FALSE))
str(brassica_SNPs)

# filter out all except major allele counts
brassica_SNPs <- brassica_SNPs[2:3388636,10:26]
colnames(brassica_SNPs) <- c("A", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "W1", "W2", "W3", 'W4', "W5", "W6", "W7", "W8")
str(brassica_SNPs)

# convert major allele counts to numeric
for (i in 1:ncol(brassica_SNPs)){
  brassica_SNPs[,i] <- with(read.table(text = brassica_SNPs[,i], sep = "/"), V1 / V2)}
str(brassica_SNPs)
brassica_SNPs_t <- as.data.frame(t(brassica_SNPs))

# making subset that I used to test
brassica_SNPs_sub <- brassica_SNPs[sample(nrow(brassica_SNPs), 100000), ]
brassica_SNPs_sub_t <- as.data.frame(t(brassica_SNPs_sub))

# run PCA
set.seed(10)
SNP_pca <- prcomp(brassica_SNPs_t, center = TRUE, scale. = TRUE)
summary(SNP_pca)

# scale pca results
x2 <- as.data.frame(scale(SNP_pca$x[,1:2]))

# reorder rows
x3 <- rbind(x2[10:17,],x2[2:9,],x2[1,])

# creating  lists for pch and col parameters 
pch.group <- c(rep(21, times=17))
col.group <- c(rep("lightskyblue", times=8), rep("tomato1", times=8), rep("gray", times=1))

# plotting Fig 1
Fig1 <- plot(x3[,1], x3[,2], xlab="PCA 1 (9.9%)", ylab="PCA 2 (9.3%)", pch=pch.group, col="black", bg=col.group, cex=1.5, las=1, asp=1)
Fig1
ggplot2::ggsave(filename = "Fig1.pdf", plot = Fig1, width = 18.4, height = 18.4, units = "cm", device='pdf', dpi=300)


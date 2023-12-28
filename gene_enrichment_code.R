### Gene enrichment


## part 3A: enrichment of SNPs from fishers exact tests (FETs) between ancestor regime and drought populations 

# fisher test for FET SNPs in exons
# 35834552 total snps with 4913424 in exons
# 913 sig snps with 147 in exons
# so 35833639 non-sig SNPs and 4913277 non-sig SNPs not in exons
fisher.test(matrix(c(35833639,4913277,913,147),byrow = TRUE, ncol = 2), alternative = "greater")
#p-value = 0.04101; differentiated SNPs between ancestors and drought descendant replicates overrepresented in exons

# fisher test for FET SNPs in Franks 2016 differentiated genes
# 35834552 total SNPs with 73528 SNPs in Franks 2016 genes
# 913 total sig snpswith 3 in Franks 2016 genes
# so 35833639 non-sig SNPs and 73525 non-sig SNPs not in Franks 2016 exons
fisher.test(matrix(c(35833639,73525,913,3),byrow = TRUE, ncol = 2), alternative = "greater")
#p-value = 0.2897; no evidence for overrepresentation of FET SNPs in Franks 2016 genes

# fisher tests for FET SNPs in flowering exons
# 35834552 total SNPs with 70352 SNPs in flowering exons
# 913 sig SNPs with 1 in flowering exon
# so 35833639 non-sig SNPs and 70351 non-sig SNPs not in flowering exons
fisher.test(matrix(c(35833639,70351,913,1),byrow = TRUE, ncol = 2), alternative = "greater")
#p-value = 0.8335; no evidence for overrepresentation of FET SNPs in flowering exons

# fisher tests for FET SNPs in stress exons
#35834552 total SNPs with 678128 SNPs in stress exons
# 913 sig SNPs with 9 in stress exons
# so 35833639 non-sig SNPs and 678119 non-sig SNPs not in stress exons
fisher.test(matrix(c(35833639,678119,913,9),byrow = TRUE, ncol = 2), alternative = "greater")
# p-value = 0.9888; no evidence for overrepresentation of FET SNPs in stress exons


## part 3B: enrichment of bayesian model snps between ancestors and descendants

# fisher text for BF > 10 SNPs in exons
# 375883 total SNPs with 49827 in exons
# 472 sig snps (BF > 10) and 57 sig snps (BF > 10) in exons
# so 375411 non-sig snps and 49772 non-sig snps not in exons

fisher.test(matrix(c(375411,49772,472,57),byrow = TRUE, ncol = 2), alternative = "greater")
#p-value = 0.7659; no evidence for overrepresentation of BF10 SNPs in  exons

# fisher text for BF > 20 SNPs in exons
# 375883 total SNPs with 49827 in exons
# 54 sig snps (BF > 20) and 11 sig snps (BF > 20) in exons
# so 375829 non-sig snps and 49818 non-sig snps not in exons

fisher.test(matrix(c(375829,49818,54,11),byrow = TRUE, ncol = 2), alternative = "greater")
#p-value = 0.1343; no evidence for overrepresentation of BF20 SNPs in exons

# fisher test fisher test BF > 10 SNPs stress genes
# 375883 total SNPs in FT genes with 6870 in stress exons
# 472 total sig SNPs (BF > 10) in stress exons and 8 in stress exons
# so 375411 non-sig snps and 6862 non-sig snps not in stress exons

fisher.test(matrix(c(375411,6862,472,8),byrow = TRUE, ncol = 2), alternative = "greater")
#p-value = 0.6312; no evidence for overrepresentation of BF10 SNPs in stress exons

# no BF > 10 SNPs identified in flowering genes to test enrichment

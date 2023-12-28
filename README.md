# README file for Johnson et al. 2023 manuscript published in Journal of Evolutionary Biology

## Article Citation:
Johnson, S. E., Tittes, S., & Franks, S. J. (2023). Rapid, nonparallel genomic evolution of Brassica rapa (field mustard) under experimental drought. Journal of Evolutionary Biology, 36(3), 550-562.

## Repository location:
DOI: 10.5061/dryad.n5tb2rbzx

## Note:
Files are names as 2022 because that is when the final analyses and submission were prepared, but the manuscript was published in 2023.

## Data files:

"Truseq3-PE-2.fa" contains adapters used to trim reads with trimmomatic software.  
"bam.txt" is a list of bam file names used in the fais_baypass.py script.  
"brassica.poolsize" is a list of the ploidy number of each population used when running the core and axillary models in BayPass.  
"brassica.covariates" is a file that encodes whether each population is ancestor or descendant for the first covariate (row 1) and, if descendant, drought or watered for the second covariate (row 2) when runing the axillary model in BayPass.  
"Johnson_et_al_2022_Evolution_individuals_data.csv" contains phenotypic data and fitness collected from the test generation and previously published in Johnson et al. 2022 in Evolution.  
"Johnson_et_al_2022_Evolution_individuals_data.xlsx" contains the same data as Johnson_et_al_2022_Evolution_individuals_data.csv in the first tab and a descriptioin of variables in the second tab.  
"Atha__v__Brap.tsv" contains a list of B. rapa genes with A. thaliana annotations added with OrthoFinder and was obtained from Prof Michael Barker at the University of Arizona.
  * Atha__v__Brap.tsv variables:
    * Orthogroup: Unique nine-digit orthogroup ID between the TAIR Arabidopsis thaliana genome (https://academic.oup.com/nar/article/40/D1/D1202/2903058) and the V1.5 Brassica rapa genome from BRAD (http://brassicadb.cn/#/SpeciesInfo/).
    * Atha: AGI locus code (AGI) consisting of the prefix "At", a one-digit chromosome number, the character "G" for gene, a unique five-digit number, and a period followed by a one digit number describing the splice variant.
    * Brap: Unique identifier code for Brassica rapa genes consisting of the prefix "Bra" and then a unique six-digit number.

## Code files:

"Johnson_et_al_2022_JEB_code.md" contains all command line code used to generate published results.  
"Fig1_code.R" contains all R code used to generate Figure 1, which plots snp frequencies of populations on pc axis.  
"Fig2_code.R" contains all R code used to generate Figure 2, which displays manhattan plots of differentiated snps identified 1) between ancestors and decendants with our Bayesian model in BayPass and 2) between drought replicates and ancestos with fisher exact tests.  
"tajimas_pi_code.R" contains all R code used to calculate mean and se nucleotide diversity (pi) from .pi files generated with PoPoolation software.  
"gene_enrichment_code.R" contains all R code used to run enrichment analysis in R from counts of significant and nonsignificant SNPs generated with bedtools intersect.  
"effective_pop_size_estimate_code.R" contains all R code used to estimate effective population size from reproduction data provided in Johnson_et_al_2022_Evolution_individuals_data.csv.  
"fais_baypass.py" produces a baypass format .txt file for downstream Baypass analysis using the reference genome and a .txt file with the names of the sorted bam files to use.  
"filter_baypass.R" is used to filter close variable loci (if within 100 bp, keeping the one with most variation) prior to running models in BayPass.  
"build_genepairs.py" contains a python script to make a table of B. rapa and A. thaliana gene pairs from Atha__v__Brap.tsv.  
"build_genepairs_tests.py" contains a python script to run unit tests.  
"Flowering_genes.R" created flowering_genes.tsv list from http://brassicadb.cn/#/FlowerGene/  
"Stress_genes.R" created stress_genes.tsv list from NCBI list of stress genes.  
"merge_stress_flowering.R" merges flowering_genes.tsv and stress_genes.tsv to create flowering_and_stress_genes.tsv.  

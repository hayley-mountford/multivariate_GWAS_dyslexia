# Script 7.1: format sumstats for SBayesRC
library(tidyverse)

# setwd

# Sumstats formats - in tab format
# * SNP: SNP name to match LD reference; A1: effect allele; A2: alternative allele; 
# * freq: frequencies of A1; b: marginal effects (reference to A1); se: sd of marginal effects; 
# * p: p value; N: per-SNP sample size;

# set sample sizes in the file as N = 1228832 for trait 1

# Format: SNP A1 A2 freq b se p N


# Import sumstats - trait 1 dyslexia with Odds ratio
dyslexia <- read_tsv("mtag_reading_GCOFF_EUR_trait_1_with_OR.txt.gz")

dyslexia_df <- select(dyslexia, c("SNP", "A1", "A2", "FRQ", "mtag_beta", "mtag_se", "mtag_pval")) 

dyslexia_df <- dyslexia_df %>% rename(p = mtag_pval, se = mtag_se, b = mtag_beta, freq = FRQ)

# Add in N
dyslexia_df$N <- 1228832

write_tsv(dyslexia_df, "mtag_reading_GCOFF_EUR_trait_1_sbayes.txt")

# Done
library(tidyverse)

# setwd

# Read MTAG trait 1 DYSLEXIA sumstats with OR in for dyslexia
sumstats <- read_tsv("mtag_reading_GCOFF_EUR_trait_1.txt.gz")

# SumHer/ LDAK requires OR for binary GWAS so generated file wit OR not Z. 
sumstats$mtag_logodds <- sumstats$mtag_z * sumstats$mtag_se 
sumstats$mtag_OR <- exp(sumstats$mtag_logodds)

range(sumstats$mtag_OR) # looks sensible 0.9705 - 1.029

write_tsv(sumstats, "mtag_reading_GCOFF_EUR_trait_1_with_OR.txt.gz")

#---------------------------------------------------------------------------------------------------#
# 4.2 FORMATTING FOR SUMHER/ LDAK - Trait 1 Dyslexia
#---------------------------------------------------------------------------------------------------#
# Sumstat formatting to run
#Predictor - the name of the predictor (ideally in the form Chr:BP, see below).
#A1 - the test allele (must be a single character, e.g., A, C, G or T).
#A2 - the other allele (must also be a single character, e.g., A, C, G or T).
#n - number of samples used when testing the predictor - need to change this to the effective pop size, not what's in this file!
#Direction - indicates whether the effect is positive or negative (with respect to the test allele). This can be an estimate of effect size (log odds).
#P - provides a p-value

#Some other things:	
#- When using SNP data, we recommend excluding predictors with ambiguous alleles (A&T or C&G), in order to avoid strand errors.
#- Only single-character alleles are allowed (usually A, C, G and T)
#- Use Chr:BP instead of rsids

df <- read_tsv("mtag_reading_GCOFF_EUR_trait_1_with_OR.txt.gz")


# Predictor chr:bp
df$Predictor <- paste(df$CHR, df$BP, sep = ":")

# Make n as N = 1228832
df$n <- 1228832

# A1 and A2 are fine

# Rename Direction and P
df <- df %>% rename(Direction = mtag_logodds, P = mtag_pval)

# Select six fields
ldak_df <- df %>% select(Predictor, A1, A2, n, Direction, P)


# Check for only allele characters - yes, only ACTG
table(ldak_df$A1)
table(ldak_df$A2)

# Check for ambiguous alleles (A&T or C&G)
dyslexia_unambig <- ldak_df %>% filter(! (A1 != "A" & A2 == "T")) %>%
  filter(! (A1 == "T" & A2 == "A")) %>%
  filter(! (A1 == "C" & A2 == "G")) %>% 
  filter(! (A1 == "G" & A2 == "C"))
#  N SNPs = 5449985 so no unambiguous SNPs :-)

# Write SumHer/LDAK formatted sumstats and zip
write_tsv(dyslexia_unambig, "mtag_reading_GCOFF_EUR_trait_1_ldak_format.txt.gz")

#---------------------------------------------------------------------------------------------------#
# DONE
#---------------------------------------------------------------------------------------------------#

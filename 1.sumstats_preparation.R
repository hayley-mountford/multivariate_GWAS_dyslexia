## R script for preparing sumstats for meta-analysis
# Load packages
library(tidyverse)

#---------------------------------------------------------------------------------------------------#
# 1.1 SUMMARY STATISTICS PREPARATION - 23andMe Dyslexia
#---------------------------------------------------------------------------------------------------#
# File is Europeans only without GC correction and filtered for imputation quality provided by 23andMe
# Read in 23andMe filtered sumstats file: dyslexia.filtered.dat 
df <- read_delim("dyslexia.filtered.dat")

# Drop chr from "scaffold" column
df$scaffold <- str_remove(df$scaffold, "chr")

# Split alleles by / and then switch order to A1/ A2
df <- separate(df, alleles, into=c("a2", "a1"), sep="/")

# Drop rows with D/I genotypes to keep only ATCG sites (MTAG requires this)
df <- df[!grepl ("D", df$a2),]

# Generate new col. called n by calculating n = im.num.0 + im.num.1 (cases + controls)
df <- mutate(df, n = im.num.0 + im.num.1)

# Used dose.b (imputed dosage estimates) as AFs
# Reorder cols and drop ones I don't want
df <- select(df, assay.name, scaffold, position, a1, a2, dose.b, effect, stderr, pvalue, n)


# Rename headers to match MTAG requirements - effect is log odds
dyslexia_23andme <- rename(df, snpid = assay.name, chr = scaffold, bpos = position, freq = dose.b, 'pval' = "pvalue") 

write_tsv(dyslexia_23andme, "dyslexia_23andme.uncorrected.txt")


# Convert to Z score - Z = effect/SE
dyslexia_23andme_zscore <- mutate(dyslexia_23andme, z = effect / stderr, .keep = "all") %>%
  select(snpid, chr, bpos, a1, a2, freq, z, pval, n)

write_tsv(dyslexia_23andme_zscore, "dyslexia_23andme_zscores_uncorrected.txt")

#---------------------------------------------------------------------------------------------------#
# 1.2 SUMMARY STATISTICS PREPARATION - GenLang Word Reading
#---------------------------------------------------------------------------------------------------#
# File is Europeans only without GC correction
# Available from https://archive.mpi.nl/mpi/islandora/object/mpi%3A1839_9058c092_111a_448c_92b7_143c434c8ea0 
# METAANALYSIS_WR_RT_EUR_combined_STERR_GCOFF_1.tbl

df <- read_tsv("../raw_data/METAANALYSIS_WR_RT_EUR_combined_STERR_GCOFF_1.tbl")

# subset the cols we need
df <- select(df, MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, 'P-value', TotalSampleSize)

# rename the cols we need
df <- rename(df, snpid = MarkerName, a1 = Allele1, a2 = Allele2, effect = Effect, stderr = StdErr, freq = Freq1, 'pval' = 'P-value', n = TotalSampleSize)


# Cross ref with chr and bpos from hg19 snp151 - just SNPs, with chr and end bp and rsid cols
# Details are in the MTAG_word_reading_details word file. 

# Sumstats have rsid but not chr or bp so annotated from dbSNP151
# Read in snp151_snps.txt which is dbsnp151 (deletions,  insertions, indels excluded)
snp151 <- read_delim ("snp151_snps.txt", delim = " ", col_names = FALSE)

# rename the headers
snp151 <- rename(snp151, chr = X1, bpos = X2, snpid = X3) 

# the left_join to annotate chr and bp by snpid
df <- left_join(df, snp151_1, by = "snpid") %>%
  filter(!is.na(chr))


# Exclude chr17_ctg5_hap1 which is included then...
# Calculate frequencies of answers grouped by chromosome to check
df <- filter(df, chr != 'chr17_ctg5_hap1')

chr_total <- df %>%
  group_by(chr) %>%
  count(chr)


# remove chr from chromosome position from chr col
df$chr <- str_remove(df$chr, "chr")

# Rename headers to match MTAG requirements - effect is beta
# sort file by chr then bpos, and rename cols
df <- arrange(df, chr, bpos) %>%
  rename(effect = Effect, stderr = StdErr, 'pval' = 'P-val')

write_tsv(df, "word_reading_EUR_genlang_GCOFF.txt")


# Convert to Z score - Z = beta/SE
word_reading_zscore <- mutate(df, z = effect / stderr, .keep = "all") %>%
  select(snpid, chr, bpos, a1, a2, freq, z, pval, n)

write_tsv(word_reading_zscore, "word_reading_EUR_genlang_GCOFF_zscores.txt")

#---------------------------------------------------------------------------------------------------#
# DONE
#---------------------------------------------------------------------------------------------------#


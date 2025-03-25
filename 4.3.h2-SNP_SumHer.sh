# Shell script for h2SNP estimation using SumHer in LDAK on linux
# Following the walk through example on https://dougspeed.com/snp-heritability/

# Installed using Conda because it was easy:
# create a new environment and install ldak6
conda create -n ldak_env -c genomedk ldak6

# load the new environment
conda activate ldak_env

# check ldak runs if you type the following
ldak6


# Download precomputed taggings computed using UK Biobank data EUR pop for all SNPs
# Recommend using the BLD-LDAK Model when estimating SNP heritability 
wget https://genetics.ghpc.au.dk/doug/bld.ldak.hapmap.gbr.tagging.gz

# And unzip
gunzip bld.ldak.hapmap.gbr.tagging.gz


# Run SumHer
# 5% sample, 5% population
ldak6 --sum-hers ./prev_permutations/mtag_trait_1_h2_sample.0.05_pop.0.05 \
		--ascertainment 0.05 \
		--prevalence 0.05 \
		--cutoff 0.01 \
		--summary mtag_reading_GCOFF_EUR_trait_1_ldak_format.txt \
		--tagfile bld.ldak.hapmap.gbr.tagging \
		--check-sums NO

# 5% sample, 7% population
ldak6 --sum-hers ./prev_permutations/mtag_trait_1_h2_sample.0.05_pop.0.07 \
		--ascertainment 0.05 \
		--prevalence 0.07 \
		--cutoff 0.01 \
		--summary mtag_reading_GCOFF_EUR_trait_1_ldak_format.txt \
		--tagfile bld.ldak.hapmap.gbr.tagging \
		--check-sums NO

# 5% sample, 10% population
ldak6 --sum-hers ./prev_permutations/mtag_trait_1_h2_sample.0.05_pop.0.1 \
		--ascertainment 0.05 \
		--prevalence 0.1 \
		--cutoff 0.01 \
		--summary mtag_reading_GCOFF_EUR_trait_1_ldak_format.txt \
		--tagfile bld.ldak.hapmap.gbr.tagging \
		--check-sums NO
		
# End
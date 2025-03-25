# Shell script for estimating h2 SNP with LDSC on linux
# LDSC was installed using conda via git clone

# git clone LDSC
git clone https://github.com/bulik/ldsc.git
cd ldsc

# Install using the Anaconda env that I set up for MTAG
conda env create --file environment.yml
source activate ldsc


# Obtained ld reference panel: /ld_ref_panel/eur_w_ld
# Obtained HapMap3 SNPs: LDSCORE_w_hm3.snplist

# Analysis used the MTAG effective GWAS sample size as fixed N
# Ignores N col, and uses Z score for consistency

#--------------------------------------------------------------#
# MTAG RESULTS - TRAIT 1 DYSLEXIA
#--------------------------------------------------------------#
# Munge sumstats
conda activate ldsc

./munge_sumstats.py --sumstats ./mtag_reading_GCOFF_EUR_trait_1.txt 
		--chunksize 50000 \
		--ignore Z_1,N,BETA \
		--N 1228832 \
		--merge-alleles LDSCORE_w_hm3.snplist \
		--out ./output/mtag_EUR_GCOFF_trait_1_munged_N_z


# 5% sample, 5% population
./ldsc.py --h2 ./output/mtag_EUR_GCOFF_trait_1_munged_N_z.sumstats.gz 
		--samp-prev 0.05 --pop-prev 0.05 \
		--ref-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--w-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--out ./output/mtag_trait_1_h2_sample.0.05_pop.0.05

# 5% sample, 7% population
./ldsc.py --h2 ./output/mtag_EUR_GCOFF_trait_1_munged_N_z.sumstats.gz \
		--samp-prev 0.05 --pop-prev 0.07 \
		--ref-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--w-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--out ./output/mtag_trait_1_h2_sample.0.05_pop.0.07


# 5% sample, 10% population
./ldsc.py --h2 ./output/mtag_EUR_GCOFF_trait_1_munged_N_z.sumstats.gz \
		--samp-prev 0.05 --pop-prev 0.1 \
		--ref-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--w-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--out ./output/mtag_trait_1_h2_sample.0.05_pop.0.1


#--------------------------------------------------------------#
# MTAG RESULTS - TRAIT 2 WORD READING
#--------------------------------------------------------------#
# Munge sumstats
conda activate ldsc

./munge_sumstats.py --sumstats ./mtag_reading_GCOFF_EUR_trait_2.txt 
		--chunksize 50000 \
		--ignore Z_1,N,BETA \
		--N 102082 \
		--merge-alleles LDSCORE_w_hm3.snplist \
		--out ./output/mtag_EUR_GCOFF_trait_2_munged_N_z
		
./ldsc.py --h2 ./output/mtag_EUR_GCOFF_trait_2_munged_N_z.sumstats.gz \
		--ref-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--w-ld-chr ld_ref_panel/eur_w_ld_chr/ \
		--out ./output/mtag_trait_2_h2

#--------------------------------------------------------------#
# DONE
#--------------------------------------------------------------#
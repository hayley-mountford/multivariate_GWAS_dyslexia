# Shell script for running partitioned heritability using LDSC

# Analysis based on: https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses

# Obtain files
cts_name=Multi_tissue_chromatin
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/Multi_tissue_chromatin_1000Gv3_ldscores.tgz
tar -xvzf ${cts_name}_1000Gv3_ldscores.tgz

# Also got files:
# https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
# https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
# https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2

tar -xvzf 1000G_Phase3_baseline_ldscores.tgz


# Munge sumstats
conda activate ldsc

../ldsc/munge_sumstats.py \
		--sumstats mtag_reading_GCOFF_EUR_trait_1_with_OR.txt \
		--merge-alleles w_hm3.snplist \
		--chunksize 50000 \
		--N 1228832 \
		--out ./out/mtag_trait_1_munged


# Split ldcts file into six using sed - runs quicker:
sed -n '1,83p' Multi_tissue_chromatin.ldcts > 1.txt
sed -n '84,167p' Multi_tissue_chromatin.ldcts > 2.txt
sed -n '168,251p' Multi_tissue_chromatin.ldcts > 3.txt
sed -n '252,335p' Multi_tissue_chromatin.ldcts > 4.txt
sed -n '336,419p' Multi_tissue_chromatin.ldcts > 5.txt
sed -n '420,498p' Multi_tissue_chromatin.ldcts > 6.txt


# Run regressions for each subset:
../ldsc/ldsc.py --h2-cts ./out/mtag_trait_1_munged.sumstats.gz \
		--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
		--out mtag_dyslexia_Multi_tissue_chromatin_1 \
		--ref-ld-chr-cts 1.txt \
		--w-ld-chr weights_hm3_no_hla/weights.
		
../ldsc/ldsc.py --h2-cts ./out/mtag_trait_1_munged.sumstats.gz \
		--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
		--out mtag_dyslexia_Multi_tissue_chromatin_2 \
		--ref-ld-chr-cts 2.txt \
		--w-ld-chr weights_hm3_no_hla/weights.
		
../ldsc/ldsc.py --h2-cts ./out/mtag_trait_1_munged.sumstats.gz \
		--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
		--out mtag_dyslexia_Multi_tissue_chromatin_3 \
		--ref-ld-chr-cts 3.txt \
		--w-ld-chr weights_hm3_no_hla/weights.
		
../ldsc/ldsc.py --h2-cts ./out/mtag_trait_1_munged.sumstats.gz \
		--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
		--out mtag_dyslexia_Multi_tissue_chromatin_4 \
		--ref-ld-chr-cts 4.txt \
		--w-ld-chr weights_hm3_no_hla/weights.
		
../ldsc/ldsc.py --h2-cts ./out/mtag_trait_1_munged.sumstats.gz \
		--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. 
		--out mtag_dyslexia_Multi_tissue_chromatin_5 \
		--ref-ld-chr-cts 5.txt \
		--w-ld-chr weights_hm3_no_hla/weights.
		
../ldsc/ldsc.py --h2-cts ./out/mtag_trait_1_munged.sumstats.gz \
		--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
		--out mtag_dyslexia_Multi_tissue_chromatin_6 \
		--ref-ld-chr-cts 6.txt \
		--w-ld-chr weights_hm3_no_hla/weights.

# END
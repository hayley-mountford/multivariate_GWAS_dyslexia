# # Polygenic Index Prediction using SBayesRC

# Used the linux version for GCTB which is the overall package:
# https://cnsgenomics.com/software/gctb/#SBayesRCTutorial

# Target data prepared as reported in Bridges et al 2022 doi:10.1017/thg.2023.2
# Requires target genotype data in plink format
# Separate phenotype file: phenotypes.final.txt
# Covars are sex, array, and 10PCs


# Download Resources
# baseline model: 2.2 https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/annot_baseline2.2.zip
# LD reference (Hapmap EUR): https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbEUR_HM3.zip

# These file links above have issues, so used the Google Drive versions instead:
# downloaded from: https://drive.google.com/drive/folders/1uxnxDjRJPzo0dTpFnERS5N2NGZX5S-sU

# unzip files using unzip command and tar -xz
tar -xf ukbEUR_HM3.zip
tar -xf annot_baseline2.2.zip


## Run SBayesRC

# 1. Tidy sumstats
./gctb_2.5.2_Linux/gctb --ldm-eigen ukbEUR_HM3 \
		--gwas-summary mtag_reading_GCOFF_EUR_trait_1_sbayes.ma \
		--impute-summary \
		--out trait1 \
		--thread 4

# Check output for errors or issues - none present


# 2. Make SNP weights for PGI
# Takes >2hours 20 minutes to run
./gctb_2.5.2_Linux/gctb --ldm-eigen ukbEUR_HM3 \
    --gwas-summary trait1.imputed.ma \
    --sbayes RC \
    --annot annot_baseline2.2.txt \
    --out trait1_weights \
    --thread 4
    

# 3. Then running the actual PGI scoring with Plink v1.90b7.2
# Ran one for each phenotype, so there are in separate bfiles files with a fam for each phenotype
plink --bfile ncds_all_age7 --score trait1_weights.snpRes 2 5 8 header sum center --out ncds_age7_trait1
plink --bfile ncds_all_age11 --score trait1_weights.snpRes 2 5 8 header sum center --out ncds_age11_trait1
plink --bfile ncds_all_age16 --score trait1_weights.snpRes 2 5 8 header sum center --out ncds_age16_trait1
plink --bfile ncds_all_age23 --score trait1_weights.snpRes 2 5 8 header sum center --out ncds_age23_trait1
plink --bfile ncds_all_age33 --score trait1_weights.snpRes 2 5 8 header sum center --out ncds_age33_trait1
plink --bfile ncds_all_all_ages --score trait1_weights.snpRes 2 5 8 header sum center --out ncds_all_ages_composite_trait1

# DONE

# # Polygenic Index Prediction using PRSice2

# PRSice for multitrait GWAS of dyslexia on NCDS reading measures
# https://github.com/choishingwan/PRSice/blob/master/docs/command_detail.md

# Target data prepared as reported in Bridges et al 2022 doi:10.1017/thg.2023.2
# Requires target genotype data in plink format
# Separate phenotype file: phenotypes.final.txt
# Covars are sex, array, and 10PCs

# age7 composite reading measures
./PRSice_linux --base mtag_reading_GCOFF_EUR_trait_1_with_OR.txt \
		--A1 A1 --A2 A2 --snp SNP --chr CHR --bp BP --stat OR --pvalue P --or \
		--target ncds_all,ncds_all.fam \
		--type bed --binary-target F \
		--geno 0.05 \
		--pheno phenotypes.final.txt --pheno-col age7 \
		--cov covariates.txt --cov-col sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
		--bar-levels 0.00000005,0.00000001,0.0000005,0.0000001,0.000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
		--fastscore --score sum --perm 10000 \
		--out prs.dyslexia.score.age7

# age11 composite reading measures
./PRSice_linux --base mtag_reading_GCOFF_EUR_trait_1_with_OR.txt \
		--A1 A1 --A2 A2 --snp SNP --chr CHR --bp BP --stat OR --pvalue P --or \
		--target ncds_all,ncds_all.fam \
		--type bed --binary-target F \
		--geno 0.05 \
		--pheno phenotypes.final.txt --pheno-col age11 \
		--cov covariates.txt --cov-col sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
		--bar-levels 0.00000005,0.00000001,0.0000005,0.0000001,0.000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
		--fastscore --score sum --perm 10000 \
		--out prs.dyslexia.score.age11

# age16 composite reading measures
./PRSice_linux --base mtag_reading_GCOFF_EUR_trait_1_with_OR.txt \
		--A1 A1 --A2 A2 --snp SNP --chr CHR --bp BP --stat OR --pvalue P --or \
		--target ncds_all,ncds_all.fam \
		--type bed --binary-target F \
		--geno 0.05 \
		--pheno phenotypes.final.txt --pheno-col age16 \
		--cov covariates.txt --cov-col sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
		--bar-levels 0.00000005,0.00000001,0.0000005,0.0000001,0.000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
		--fastscore --score sum --perm 10000 \
		--out prs.dyslexia.score.age16

# age23 BINARY self-report of reading difficulties
./PRSice_linux --base mtag_reading_GCOFF_EUR_trait_1_with_OR.txt \
		--A1 A1 --A2 A2 --snp SNP --chr CHR --bp BP --stat OR --pvalue P --or \
		--target ncds_all,ncds_all.fam \
		--type bed --binary-target T \
		--geno 0.05 \
		--pheno phenotypes.final.txt --pheno-col age23 \
		--cov covariates.txt --cov-col sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
		--bar-levels 0.00000005,0.00000001,0.0000005,0.0000001,0.000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
		--fastscore --score sum --perm 10000 \
		--out prs.dyslexia.score.age23

# age33 BINARY self-report of reading difficulties
./PRSice_linux --base mtag_reading_GCOFF_EUR_trait_1_with_OR.txt \
		--A1 A1 --A2 A2 --snp SNP --chr CHR --bp BP --stat OR --pvalue P --or \
		--target ncds_all,ncds_all.fam \
		--type bed --binary-target T \
		--geno 0.05 \
		--pheno phenotypes.final.txt --pheno-col age33 \
		--cov covariates.txt --cov-col sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
		--bar-levels 0.00000005,0.00000001,0.0000005,0.0000001,0.000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
		--fastscore --score sum --perm 10000 \
		--out prs.dyslexia.score.age33

# all - composite of all ages of reading measures and self-report of reading difficulties
./PRSice_linux --base mtag_reading_GCOFF_EUR_trait_1_with_OR.txt \
		--A1 A1 --A2 A2 --snp SNP --chr CHR --bp BP --stat OR --pvalue P --or \
		--target ncds_all,ncds_all.fam \
		--type bed --binary-target F \
		--geno 0.05 \
		--pheno phenotypes.final.txt --pheno-col all \
		--cov covariates.txt --cov-col sex,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
		--bar-levels 0.00000005,0.00000001,0.0000005,0.0000001,0.000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
		--fastscore --score sum --perm 10000 \
		--out prs.dyslexia.score.all

# DONE
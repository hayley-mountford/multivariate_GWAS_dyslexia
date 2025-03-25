library(tidyverse)
library(viridis)

# setwd

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 1 - data prep
#---------------------------------------------------------------------------------------------------#
# Prep - read in covars file
covars <- read_table("covariates.txt") %>%
  select(fid, sex, array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
  rename(FID = fid)

# Read in PRSs for ages 7, 11, 17 and all (reading composite) - no p filtering
age7_prs <- read_table("ncds_age7_trait1.P1.profile", na = "-9")
age11_prs <- read_table("ncds_age11_trait1.P1.profile", na = "-9")
age16_prs <- read_table("ncds_age16_trait1.P1.profile", na = "-9")
comp_prs <- read_table("ncds_all_ages_composite_trait1.P1.profile", na = "-9")

# Update prs files to make FID
age7_prs$FID <- paste0(age7_prs$FID, "_", age7_prs$IID)
age7_prs <- age7_prs %>% select(FID, PHENO, SCORESUM)

age11_prs$FID <- paste0(age11_prs$FID, "_", age11_prs$IID)
age11_prs <- age11_prs %>% select(FID, PHENO, SCORESUM)

age16_prs$FID <- paste0(age16_prs$FID, "_", age16_prs$IID)
age16_prs <- age16_prs %>% select(FID, PHENO, SCORESUM)

comp_prs$FID <- paste0(comp_prs$FID, "_", comp_prs$IID)
comp_prs <- comp_prs %>% select(FID, PHENO, SCORESUM)

# Make a file for each binary phen with PRSs, phenotype and covars
age7_prs <- left_join(age7_prs, covars, by = 'FID')
age11_prs <- left_join(age11_prs, covars, by = 'FID')
age16_prs <- left_join(age16_prs, covars, by = 'FID')
comp_prs <- left_join(comp_prs, covars, by = 'FID')


#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 2 - calculate mean and SDs
#---------------------------------------------------------------------------------------------------#
## Check current distribs of PGIs - these are scaled
# Descriptive stats
describe(age7_prs$SCORESUM) 
## Age 7 - mean = 0, SD = 0.11

# Descriptive stats
describe(age11_prs$SCORESUM)
## Age 11 - mean = 0, SD = 0.11

# Descriptive stats
describe(age16_prs$SCORESUM)
## Age 16 - mean = 0, SD = 0.11

# Descriptive stats
describe(comp_prs$SCORESUM)
## Reading composite - mean = 0, SD = 0.11

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 3 - rescale to mean and SD of controls
#---------------------------------------------------------------------------------------------------#
# To rescale polygenic index scores (PGIs) to control the mean and standard deviation, 
# you calculate the z-score for each individual's PRS by subtracting the population mean 
# from their score and then dividing by the population standard deviation;

# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
age7_prs$PRS_scaled <- (age7_prs$SCORESUM - 0) / 0.11

# Quick histo to check distrib:
ggplot(age7_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, fairly normal :-) scale is mean is zero

# save the plot
ggsave("SBayesRC_age7_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Z-score = (Individual PRS - Population Mean) / Population Standard Deviation
age11_prs$PRS_scaled <- (age11_prs$SCORESUM - 0) / 0.11

# Quick histo to check distrib:
ggplot(age11_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, fairly normal :-) scale is mean is zero

# save the plot
ggsave("SBayesRC_age11_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
age16_prs$PRS_scaled <- (age16_prs$SCORESUM - 0) / 0.11

# Quick histo to check distrib:
ggplot(age16_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, fairly normal :-) scale is mean is zero

# save the plot
ggsave("SBayesRC_age16_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
comp_prs$PRS_scaled <- (comp_prs$SCORESUM - 0) / 0.11

# Quick histo to check distrib:
ggplot(comp_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, fairly normal :-) scale is mean is zero

# save the plot
ggsave("SBayesRC_composite_dyslexia_scaled_PRS.png", height = 7, width = 7)

# Done - next to regressions

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 4A - AGE 7 - NULL Regressions with just covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 5712
age7_for_regression <- age7_prs %>% filter(!is.na(PHENO))


## LM for age 7 - no PRS only covars
age7_null_model <- lm(PHENO ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age7_for_regression)
summary(age7_null_model)

age7.null <- summary(age7_null_model)$r.squared %>% as.data.frame()
age7.null$Thresh <- "prs"
age7.null <- age7.null %>% rename(null.r2 = ".")

## LM for age 7 - all covariates
age7_full_model <- lm(PHENO ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age7_for_regression)
summary(age7_full_model)
age7.full <- "prs" %>% as.data.frame()
age7.full$model.r2 <- summary(age7_full_model)$r.squared
age7.full$Estimate <- summary(age7_full_model)$coefficients["PRS_scaled", "Estimate"]
age7.full$se <- summary(age7_full_model)$coefficients["PRS_scaled", "Std. Error"]
age7.full$Pval <- summary(age7_full_model)$coefficients["PRS_scaled", "Pr(>|t|)"]
age7.full <- age7.full %>% rename(Thresh = ".")

# Join null and full models
age7_results <- left_join(age7.null, age7.full, by = "Thresh")
age7_results <- age7_results %>% select(Thresh, null.r2, model.r2, Estimate, se, Pval)


## Calculate R2 for the model
# R2 of PRS is calculated as the model R2 minus the null R2
age7_results$prs.r2 <- age7_results$model.r2 - age7_results$null.r2

# Write files
write_csv(age7_results, "SBayesRC_age7_dyslexia_PRS.csv")

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 4B - AGE 11 - Regressions with and without covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 5528
age11_for_regression <- age11_prs %>% filter(!is.na(PHENO))

## LM for age 11 - no PRS only covars
age11_null_model <- lm(PHENO ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age11_for_regression)
summary(age11_null_model)

age11.null <- summary(age11_null_model)$r.squared %>% as.data.frame()
age11.null$Thresh <- "prs"
age11.null <- age11.null %>% rename(null.r2 = ".")

## LM for age 11 - all covariates
age11_full_model <- lm(PHENO ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age11_for_regression)
summary(age11_full_model)
age11.full <- "prs" %>% as.data.frame()
age11.full$model.r2 <- summary(age11_full_model)$r.squared
age11.full$Estimate <- summary(age11_full_model)$coefficients["PRS_scaled", "Estimate"]
age11.full$se <- summary(age11_full_model)$coefficients["PRS_scaled", "Std. Error"]
age11.full$Pval <- summary(age11_full_model)$coefficients["PRS_scaled", "Pr(>|t|)"]
age11.full <- age11.full %>% rename(Thresh = ".")

# Join null and full models
age11_results <- left_join(age11.null, age11.full, by = "Thresh")
age11_results <- age11_results %>% select(Thresh, null.r2, model.r2, Estimate, se, Pval)

## Calculate R2 for the model
# R2 of PRS is calculated as the model R2 minus the null R2
age11_results$prs.r2 <- age11_results$model.r2 - age11_results$null.r2

# Write files
write_csv(age11_results, "SBayes_age11_dyslexia_PRS.csv")

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 4C - AGE 16 - Regressions with and without covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 4809
age16_for_regression <- age16_prs %>% filter(!is.na(PHENO))

## LM for age 16 - no PRS only covars
age16_null_model <- lm(PHENO ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age16_for_regression)
summary(age16_null_model)

age16.null <- summary(age16_null_model)$r.squared %>% as.data.frame()
age16.null$Thresh <- "prs"
age16.null <- age16.null %>% rename(null.r2 = ".")

## LM for age 16 - all covariates
age16_full_model <- lm(PHENO ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age16_for_regression)
summary(age16_full_model)
age16.full <- "prs" %>% as.data.frame()
age16.full$model.r2 <- summary(age16_full_model)$r.squared
age16.full$Estimate <- summary(age16_full_model)$coefficients["PRS_scaled", "Estimate"]
age16.full$se <- summary(age16_full_model)$coefficients["PRS_scaled", "Std. Error"]
age16.full$Pval <- summary(age16_full_model)$coefficients["PRS_scaled", "Pr(>|t|)"]
age16.full <- age16.full %>% rename(Thresh = ".")

# Join null and full models
age16_results <- left_join(age16.null, age16.full, by = "Thresh")
age16_results <- age16_results %>% select(Thresh, null.r2, model.r2, Estimate, se, Pval)

## Calculate R2 for the model
# R2 of PRS is calculated as the model R2 minus the null R2
age16_results$prs.r2 <- age16_results$model.r2 - age16_results$null.r2

# Write files
write_csv(age16_results, "SBayesRC_age16_dyslexia_PRS.csv")

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 4D - READING ALL AGES COMPOSITE - Regressions with and without covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 3089
comp_for_regression <- comp_prs %>% filter(!is.na(PHENO))

## LM for all ages composite - no PRS only covars
comp_null_model <- lm(PHENO ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = comp_for_regression)
summary(comp_null_model)

comp.null <- summary(comp_null_model)$r.squared %>% as.data.frame()
comp.null$Thresh <- "prs"
comp.null <- comp.null %>% rename(null.r2 = ".")

## LM for all ages composite - all covariates
comp_full_model <- lm(PHENO ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = comp_for_regression)
summary(comp_full_model)
comp.full <- "prs" %>% as.data.frame()
comp.full$model.r2 <- summary(comp_full_model)$r.squared
comp.full$Estimate <- summary(comp_full_model)$coefficients["PRS_scaled", "Estimate"]
comp.full$se <- summary(comp_full_model)$coefficients["PRS_scaled", "Std. Error"]
comp.full$Pval <- summary(comp_full_model)$coefficients["PRS_scaled", "Pr(>|t|)"]
comp.full <- comp.full %>% rename(Thresh = ".")

# Join null and full models
comp_results <- left_join(comp.null, comp.full, by = "Thresh")
comp_results <- comp_results %>% select(Thresh, null.r2, model.r2, Estimate, se, Pval)

## Calculate R2 for the model
# R2 of PRS is calculated as the model R2 minus the null R2
comp_results$prs.r2 <- comp_results$model.r2 - comp_results$null.r2

# Write files
write_csv(comp_results, "SBayesRC_comp_dyslexia_PRS.csv")

##--------------------------------------------------------------------------------##
## SBayesRC quantitatives are done
##--------------------------------------------------------------------------------##
library(tidyverse)
library(psych)

# setwd

#---------------------------------------------------------------------------------------------------#
# 6.3 ## Part 1 - data prep
#---------------------------------------------------------------------------------------------------#
# Prep - read in phenotype and covars file
phens <- read_table("phenotypes.final.txt") %>%
  select(fid, age7, age11, age16, all) %>%
  rename(FID = fid)

covars <- read_table("covariates.txt") %>%
  select(fid, sex, array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
  rename(FID = fid)


# Read in scores for ages 7, 11, 17 and all (reading composite)
age7_prs <- read_table("prs.dyslexia.score.age7.best")
age11_prs <- read_table("prs.dyslexia.score.age11.best")
age16_prs <- read_table("prs.dyslexia.score.age16.best")
comp_prs <- read_table("prs.dyslexia.score.all.best")

# Make a file for each binary phen with PGI, phenotype and covars
age7_prs <- left_join(age7_prs, phens, by = 'FID')
age7_prs <- left_join(age7_prs, covars, by = 'FID') %>%
  select(FID, PRS, age7, sex, array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

age11_prs <- left_join(age11_prs, phens, by = 'FID')
age11_prs <- left_join(age11_prs, covars, by = 'FID') %>%
  select(FID, PRS, age11, sex, array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

age16_prs <- left_join(age16_prs, phens, by = 'FID')
age16_prs <- left_join(age16_prs, covars, by = 'FID') %>%
  select(FID, PRS, age16, sex, array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

comp_prs <- left_join(comp_prs, phens, by = 'FID')
comp_prs <- left_join(comp_prs, covars, by = 'FID') %>%
  select(FID, PRS, all, sex, array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

#---------------------------------------------------------------------------------------------------#
# 6.3 ## Part 2 - calculate mean and SDs
#---------------------------------------------------------------------------------------------------#
## Check current distribs of PGIs - not scaled
# Descriptive stats
describe(age7_prs$PRS) 
## Age 7 - mean = -1.32, SD = 0.39

# Descriptive stats
describe(age11_prs$PRS)
## Age 11 - mean = -1.68, SD = 0.45

# Descriptive stats
describe(age16_prs$PRS)
## Age 16 - mean = -1.32, SD = 0.39

# Descriptive stats
describe(comp_prs$PRS)
## Reading composite - mean = -1.32, SD = 0.39

#---------------------------------------------------------------------------------------------------#
# 6.3 ## Part 3 - Rescale to mean and SD of controls
#---------------------------------------------------------------------------------------------------#
# To rescale polygenic index (PGI) to control the mean and standard deviation, 
# you calculate the z-score for each individual's PGI by subtracting the population mean 
# from their score and then dividing by the population standard deviation;

# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
age7_prs$PRS_scaled <- (age7_prs$PRS - -1.32) / 0.39

# Quick histo to check distrib:
ggplot(age7_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, normal and scale is mean is zero

# save the plot
ggsave("PRSice2_age7_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
age11_prs$PRS_scaled <- (age11_prs$PRS - -1.68) / 0.45

# Quick histo to check distrib:
ggplot(age11_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, fairly normal :-) scale is mean is zero

# save the plot
ggsave("PRSice2_age11_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
age16_prs$PRS_scaled <- (age16_prs$PRS - -1.32) / 0.39

# Quick histo to check distrib:
ggplot(age16_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, fairly normal :-) scale is mean is zero

# save the plot
ggsave("PRSice2_age16_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
comp_prs$PRS_scaled <- (comp_prs$PRS - -1.32) / 0.39

# Quick histo to check distrib:
ggplot(comp_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic() # Yes, fairly normal :-) scale is mean is zero

# save the plot
ggsave("PRSice2_composite_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Done - next to regressions

#---------------------------------------------------------------------------------------------------#
# 6.3 ## Part 4A - AGE 7 - NULL Regressions with just covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 5712
age7_for_regression <- age7_prs %>% filter(!is.na(age7))


## LM for age 7 - no PRS only covars
age7_null_model <- lm(age7 ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age7_for_regression)
summary(age7_null_model)

age7.null <- summary(age7_null_model)$r.squared %>% as.data.frame()
age7.null$Thresh <- "prs.0.005"
age7.null <- age7.null %>% rename(null.r2 = ".")


## LM for age 7 - all covariates
age7_full_model <- lm(age7 ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age7_for_regression)
summary(age7_full_model)
age7.full <- "prs.0.005" %>% as.data.frame()
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
write_csv(age7_results, "PRSice2_age7_dyslexia_PRS.csv")

#---------------------------------------------------------------------------------------------------#
# 6.3 ## Part 4B - AGE 11 - Regressions with and without covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 5528
age11_for_regression <- age11_prs %>% filter(!is.na(age11))

## LM for age 23 - no PRS only covars
age11_null_model <- lm(age11 ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age11_for_regression)
summary(age11_null_model)

age11.null <- summary(age11_null_model)$r.squared %>% as.data.frame()
age11.null$Thresh <- "prs.0.01"
age11.null <- age11.null %>% rename(null.r2 = ".")

## LM for age 7 - all covariates
age11_full_model <- lm(age11 ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age11_for_regression)
summary(age11_full_model)
age11.full <- "prs.0.01" %>% as.data.frame()
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
write_csv(age11_results, "PRSice2_age11_dyslexia_PRS.csv")

#---------------------------------------------------------------------------------------------------#
# 6.3 ## Part 4C - AGE 16 - Regressions with and without covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 5528
age16_for_regression <- age16_prs %>% filter(!is.na(age16))

## LM for age 23 - no PRS only covars
age16_null_model <- lm(age16 ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age16_for_regression)
summary(age16_null_model)

age16.null <- summary(age16_null_model)$r.squared %>% as.data.frame()
age16.null$Thresh <- "prs.0.005"
age16.null <- age16.null %>% rename(null.r2 = ".")

## LM for age 7 - all covariates
age16_full_model <- lm(age16 ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = age16_for_regression)
summary(age16_full_model)
age16.full <- "prs.0.005" %>% as.data.frame()
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
write_csv(age16_results, "PRSice2_age16_dyslexia_PRS.csv")

#---------------------------------------------------------------------------------------------------#
# 6.3 ## Part 4D - READING COMPOSITE - Regressions with and without covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 3089
comp_for_regression <- comp_prs %>% filter(!is.na(all))

## LM for age 23 - no PRS only covars
comp_null_model <- lm(all ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = comp_for_regression)
summary(comp_null_model)

comp.null <- summary(comp_null_model)$r.squared %>% as.data.frame()
comp.null$Thresh <- "prs.0.005"
comp.null <- comp.null %>% rename(null.r2 = ".")

## LM for age 7 - all covariates
comp_full_model <- lm(all ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = comp_for_regression)
summary(comp_full_model)
comp.full <- "prs.0.005" %>% as.data.frame()
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
write_csv(comp_results, "PRSice2_comp_dyslexia_PRS.csv")

##--------------------------------------------------------------------------------##
## PRSice2 quantitatives are done
##--------------------------------------------------------------------------------##

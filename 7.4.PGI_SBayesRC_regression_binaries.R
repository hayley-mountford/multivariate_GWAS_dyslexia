library(tidyverse)
library(viridis)
library(ISLR) # for calculating r2 from glm

# setwd


#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 1 - data prep
#---------------------------------------------------------------------------------------------------#
# Prep - read in covars file
covars <- read_table("covariates.txt") %>%
  select(fid, sex, array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
  rename(FID = fid)


# Read in PGIs for age 23 and 33 - no p filtering
age23_prs <- read_table("ncds_age23_trait1.P1.profile", na = "-9")
age33_prs <- read_table("ncds_age33_trait1.P1.profile", na = "-9")

# Update PGI files to make FID
age23_prs$FID <- paste0(age23_prs$FID, "_", age23_prs$IID)
age23_prs <- age23_prs %>% select(FID, PHENO, SCORESUM)

age33_prs$FID <- paste0(age33_prs$FID, "_", age33_prs$IID)
age33_prs <- age33_prs %>% select(FID, PHENO, SCORESUM)


# Make a file for each binary phen with PPGI, phenotype and covars
age23_prs <- left_join(age23_prs, covars, by = 'FID')
age33_prs <- left_join(age33_prs, covars, by = 'FID')


#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 2 - check distributions, mean and SD
#---------------------------------------------------------------------------------------------------#
## Check current distribs of PGIs - scaled
# Descriptive stats - Mean = 0, SD = 0.11. 
describe(age23_prs$SCORESUM)
table(age23_prs$PHENO)

# Check just controls descriptive stats - N = 5388
age23_prs_controls <- age23_prs %>% filter(PHENO == 1)
describe(age23_prs_controls$SCORESUM) 

## Age 23 controls - mean = 0, SD = 0.11


# Descriptive stats - Mean = 0, SD = 0.11
describe(age33_prs$SCORESUM)
table(age33_prs$PHENO) # 5700 in cohort

# Check just controls descriptive stats - N = 5388
age33_prs_controls <- age33_prs %>% filter(PHENO == 1)
describe(age33_prs_controls$SCORESUM) 

## Age 33 controls - mean = 0, SD = 0.11

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 3 - rescale to mean and SD of controls
#---------------------------------------------------------------------------------------------------#
# To rescale polygenic index (PGIs) to control the mean and standard deviation, 
# you calculate the z-score for each individual's PGI by subtracting the population mean 
# from their score and then dividing by the population standard deviation;

# Z-score = (Individual PGI - Population Mean) / Population Standard Deviation
age23_prs$PRS_scaled <- (age23_prs$SCORESUM - 0) / 0.11

# Quick histo to check distrib:
ggplot(age23_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) # Yes, fairly normal :-) scale is mean is zero

# Box plot of PRS for dyslexia on reading difficulties at age 23 by case/ control
ggplot(age23_prs, aes(x = PHENO, y = PRS_scaled, color = PHENO)) + 
  scale_x_discrete(limits=c("1", "2"), labels = c('No','Yes')) +
  scale_y_continuous(breaks = seq(-5, 5, by = 1)) +
  labs(x="Difficulties with reading at age 23", y="Dyslexia PRS") +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend = FALSE) +
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2, color="black" ) +
  theme_classic()

# save the plot
ggsave("SBayesRC_age23_dyslexia_scaled_PRS.png", height = 7, width = 7)


# Now for age 33
age33_prs$PRS_scaled <- (age33_prs$SCORESUM - 0) / 0.11

# Quick histo to check distrib:
ggplot(age33_prs, aes(x = PRS_scaled)) +
  geom_histogram(binwidth = 0.05) # Yes, fairly normal :-) scale is mean is zero

# Box plot of PRS for dyslexia on reading difficulties at age 33 by case/ control
ggplot(age33_prs, aes(x = PHENO, y = PRS_scaled, color = PHENO)) + 
  scale_x_discrete(limits=c("1", "2"), labels = c('No','Yes')) +
  scale_y_continuous(breaks = seq(-4, 4, by = 1)) +
  labs(x="Difficulties with reading at age 33", y="Dyslexia PRS") +
  geom_jitter(shape = 16, position = position_jitter(0.2), show.legend = FALSE) +
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2, color="black" ) +
  theme_classic()

# save the plot
ggsave("SBayes_age33_dyslexia_scaled_PRS.png", height = 7, width = 7)

# Done - next to regressions

#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 4A - AGE 23 - Regressions with covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 5555
age23_for_regression <- age23_prs %>% filter(!is.na(PHENO))

# Recode 0 as control and 1 as case
age23_for_regression$PHENO[age23_for_regression$PHENO == 1] <- 0
age23_for_regression$PHENO[age23_for_regression$PHENO == 2] <- 1


## GLM for age 23 - no PRS only covars
age23_null_model <- glm(PHENO ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = age23_for_regression)
summary(age23_null_model)

# calculate McFadden's R-squared for model using ISLR package
age23_null <- with(summary(age23_null_model), 1 - deviance/null.deviance) %>% as.data.frame()
age23_null$Thresh <- "prs"
age23_null <- age23_null %>% rename(null.r2 = ".")


## GLM for age 23 - all covariates
age23_full <- glm(PHENO ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = age23_for_regression)
summary(age23_full)

results.age23.full <- "prs" %>% as.data.frame()
results.age23.full$model.r2 <- with(summary(age23_full), 1 - deviance/null.deviance)
results.age23.full$Estimate <- summary(age23_full)$coefficients["PRS_scaled", "Estimate"]
results.age23.full$se <- summary(age23_full)$coefficients["PRS_scaled", "Std. Error"]
results.age23.full$Pval <- summary(age23_full)$coefficients["PRS_scaled", "Pr(>|z|)"]
results.age23.full <- results.age23.full %>% rename(Thresh = ".")

# Join null and full models
age23_results <- left_join(age23_null, results.age23.full, by = "Thresh")
age23_results <- age23_results %>% select(Thresh, null.r2, model.r2, Estimate, se, Pval)


## Calculate R2 for the model
# R2 of PRS is calculated as the model R2 minus the null R2
age23_results$prs.r2 <- age23_results$model.r2 - age23_results$null.r2
age23_results$Pval <- signif(age23_results$Pval, digits=3)

## Calculate Odds ratio and CIs
# Odds Ratio (OR): OR = exp(Effect)
age23_results$OR <- exp(age23_results$Estimate)

# Confidence Interval (CI): [exp(Effect - 1.96 * SE), exp(Effect + 1.96 * SE)] 
age23_results$CI_lower <- exp(age23_results$Estimate - 1.96 * results.age23.full$se)
age23_results$CI_upper <- exp(age23_results$Estimate + 1.96 * results.age23.full$se)


# Write files
write_csv(age23_results, "SBayesRC_age23_dyslexia_PRS.csv")


#---------------------------------------------------------------------------------------------------#
# 7.3 ## Part 4B - AGE 33 - Regressions with and without covars
#---------------------------------------------------------------------------------------------------#
# Select just people with phenotypes N = 5555
age33_for_regression <- age33_prs %>% filter(!is.na(PHENO))

# Recode 0 as control and 1 as case
age33_for_regression$PHENO[age33_for_regression$PHENO == 1] <- 0
age33_for_regression$PHENO[age33_for_regression$PHENO == 2] <- 1


## GLM for age 23 - no PRS, only covars
age33_null_model <- glm(PHENO ~ sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = age33_for_regression)
summary(age33_null_model)

# calculate McFadden's R-squared for model using ISLR package
age33_null <- with(summary(age33_null_model), 1 - deviance/null.deviance) %>% as.data.frame()
age33_null$Thresh <- "prs"
age33_null <- age33_null %>% rename(null.r2 = ".")


## GLM for age 23 - all covariates
age33_full <- glm(PHENO ~ PRS_scaled + sex + array + PC1 + PC2 + PC3 + PC4  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial, data = age33_for_regression)
summary(age33_full)

results.age33.full <- "prs" %>% as.data.frame()
results.age33.full$model.r2 <- with(summary(age33_full), 1 - deviance/null.deviance)
results.age33.full$Estimate <- summary(age33_full)$coefficients["PRS_scaled", "Estimate"]
results.age33.full$se <- summary(age33_full)$coefficients["PRS_scaled", "Std. Error"]
results.age33.full$Pval <- summary(age33_full)$coefficients["PRS_scaled", "Pr(>|z|)"]
results.age33.full <- results.age33.full %>% rename(Thresh = ".")

# Join null and full models
age33_results <- left_join(age33_null, results.age33.full, by = "Thresh")
age33_results <- age33_results %>% select(Thresh, null.r2, model.r2, Estimate, se, Pval)


## Calculate R2 for the model
# R2 of PRS is calculated as the model R2 minus the null R2
age33_results$prs.r2 <- age33_results$model.r2 - age33_results$null.r2
age33_results$Pval <- signif(age33_results$Pval, digits=3)

## Calculate Odds ratio and CIs
# Odds Ratio (OR): OR = exp(Effect)
age33_results$OR <- exp(age33_results$Estimate)

# Confidence Interval (CI): [exp(Effect - 1.96 * SE), exp(Effect + 1.96 * SE)] 
age33_results$CI_lower <- exp(age33_results$Estimate - 1.96 * results.age33.full$se)
age33_results$CI_upper <- exp(age33_results$Estimate + 1.96 * results.age33.full$se)


# Write files
write_csv(age33_results, "SBayes_age33_dyslexia_PRS.csv")

##--------------------------------------------------------------------------------##
## PRSice2 binaries are done
##--------------------------------------------------------------------------------##
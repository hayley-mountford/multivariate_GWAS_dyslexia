library(tidyverse)
library(shadowtext)

# setwd

# Read in MTAG sumstats
dyslexia_raw <- read_tsv("mtag_reading_GCOFF_EUR_trait_1.txt")
reading_raw <- read_tsv("mtag_reading_GCOFF_EUR_trait_2.txt")

#---------------------------------------------------------------------------------------------------#
# 3.1 PLOT MANHATTAN - Trait 1 Dyslexia
#---------------------------------------------------------------------------------------------------#
# Based on https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
# Thinning settings for final plot 5e-2 (sig) and 0.2 (slice)
# Thinning settings for testing 1e-5 (sig) and 0.001 (slice)

# Datapoint thinning (because 5.5M SNPs is too many)
sig_data <- dyslexia_raw %>% 
  subset(mtag_pval < 5e-2)

notsig_data <- dyslexia_raw %>% 
  subset(mtag_pval >= 5e-2) %>%
  group_by(CHR) %>% 
  slice_sample(prop = 0.2)

thinned_data <- bind_rows(sig_data, notsig_data) 


# Calculates cumulative position for the x axis positions, and makes a new variable (bp_cum) with this position
data_cum <- thinned_data %>%
  mutate(CHR = factor(CHR, levels = c(as.character(1:22), "X"))) %>%
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHR, bp_add)

gwas_data <- thinned_data %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

gene_names <- read_csv("gene_labels_dyslexia.csv")

gwas_data <- left_join(gwas_data, gene_names, by = "SNP", keep = FALSE)


# Set thresholds 
axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(mtag_pval == min(mtag_pval)) %>% 
  mutate(ylim = abs(floor(log10(mtag_pval))) + 2) %>% 
  pull(ylim)

sig <- 5e-8


# Import GenLang word reading GWAS
genlang <- read_tsv("word_reading_EUR_genlang_GCOFF.txt")

# filter only sig snps - NONE ARE SIG
sig_data_genlang <- genlang %>% 
  subset(pval < 5e-8) %>%
  rename(CHR = chr, BP = bpos, mtag_pval = pval, SNP = snpid)

# So instead - filter in rs11208009
sig_data_genlang <- genlang %>%
  filter(snpid == "rs11208009") %>%
  rename(CHR = chr, BP = bpos, mtag_pval = pval, SNP = snpid)

  
# Make cumulative position for GenLang SNP
gwas_data_gl <- sig_data_genlang %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)


# Import dyslexia GWAS
dyslexia <- read_tsv("dyslexia_23andme1.17623GCcorrected.txt")

# filter only sig snps
sig_data_23andme <- dyslexia %>% 
  subset(pval < 5e-8) %>%
  rename(CHR = chr, BP = bpos, mtag_pval = pval, SNP = snpid)

# add bps
gwas_data_23andme <- sig_data_23andme %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)


# Set col for colour for each dataset
just_sig_23andme <- select(gwas_data_23andme, "SNP") %>%
  mutate(colour = "p")

just_sig_genlang <- select(gwas_data_gl, "SNP") %>%
  mutate(colour = "r")
  
just_sig_snps <- rbind(just_sig_23andme, just_sig_genlang)


# Join info from GenLang and dyslexia to MTAG sumstats
gwas_data_2 <- left_join(gwas_data, just_sig_snps, by = "SNP")


# Then plot
mtag_dyslexia_plot <- ggplot(gwas_data_2, aes(x = bp_cum, y = -log10(mtag_pval), color = as_factor(CHR))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75, show.legend = FALSE) +
  geom_point(data = gwas_data_2 %>% filter(colour == "p"), color = "purple") +
  geom_point(data = gwas_data_2 %>% filter(colour =="r"), color = "red") +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log10(P-value)") + 
  theme_minimal() +
  theme(
    legend.position="none", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  geom_shadowtext(data = subset(gwas_data_2, !is.na(gene)), 
    aes(label = gene), hjust = "left", angle = 45, color = "black", bg.colour = "white")

ggsave("manhattan_plot_dyslexia.png", mtag_dyslexia_plot, bg="white", height = 5, width = 8, scale = 2.2)


#---------------------------------------------------------------------------------------------------#
# 3.2 PLOT MANHATTAN - Trait 2 Word Reading
#---------------------------------------------------------------------------------------------------#
# Datapoint thinning (because 5.5M SNPs is too many)
sig_data <- reading_raw %>% 
  subset(mtag_pval < 1e-3)

notsig_data <- reading_raw %>% 
  subset(mtag_pval >= 1e-3) %>%
  group_by(CHR) %>% 
  slice_sample(prop = 0.1)

thinned_data <- bind_rows(sig_data, notsig_data) 


# Calculates cumulative position for the x axis positions, and makes a new variable (bp_cum) with this position
data_cum <- thinned_data %>%
  mutate(CHR = factor(CHR, levels = c(as.character(1:22), "X"))) %>%
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHR, bp_add)

gwas_data <- thinned_data %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

# Set thresholds 
axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(mtag_pval == min(mtag_pval)) %>% 
  mutate(ylim = abs(floor(log10(mtag_pval))) + 2) %>% 
  pull(ylim)

sig <- 5e-8

# Then plot # aes size = -log10(mtag_pval)
word_reading_plot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(mtag_pval), color = as_factor(CHR))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75, show.legend = FALSE) +
  geom_point(data = gwas_data_23andme, color = "purple") +
  geom_point(data = gwas_data_gl, color = "red") +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$CHR)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log10(P-value)") + 
  theme_minimal() +
  theme(
    legend.position="none", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave("manhattan_plot_word_reading.png", mtag_word_reading_plot, bg="white", height = 5, width = 8, scale = 2.2)

#---------------------------------------------------------------------------------------------------#
# DONE
#---------------------------------------------------------------------------------------------------#

# Script for partitioned heritability results from LDSC
library(tidyverse)
library(viridis)

# setwd

#---------------------------------------------------------------------------------------------------#
# 5.2 Import partitioned heritability results - Trait 1 Dyslexia
#---------------------------------------------------------------------------------------------------#
# Import txt files and append to new df
df_1 <- read_tsv("mtag_dyslexia_Multi_tissue_chromatin_1.cell_type_results.txt") 
df_2 <- read_tsv("mtag_dyslexia_Multi_tissue_chromatin_2.cell_type_results.txt")
df_3 <- read_tsv("mtag_dyslexia_Multi_tissue_chromatin_3.cell_type_results.txt")
df_4 <- read_tsv("mtag_dyslexia_Multi_tissue_chromatin_4.cell_type_results.txt")
df_5 <- read_tsv("mtag_dyslexia_Multi_tissue_chromatin_5.cell_type_results.txt")
df_6 <- read_tsv("mtag_dyslexia_Multi_tissue_chromatin_6.cell_type_results.txt")

df <- bind_rows(df_1, df_2, df_3, df_4, df_5, df_6)

#---------------------------------------------------------------------------------------------------#
# 5.2 Group and write full results table - Trait 1 Dyslexia
#---------------------------------------------------------------------------------------------------#
# Group tissues
df$Tissue <- "Other"
df$Tissue[grep("Brain", df$Name, ignore.case = TRUE)] <- "Brain"
df$Tissue[grep("Muscle", df$Name, ignore.case = TRUE)] <- "Muscle"
df$Tissue[grep("Blood", df$Name, ignore.case = TRUE)] <- "lmmune or blood"
df$Tissue[grep("Thymus", df$Name, ignore.case = TRUE)] <- "lmmune or blood"
df$Tissue[grep("Spleen", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("skin", df$Name, ignore.case = TRUE)] <- "Skin"
df$Tissue[grep("neur", df$Name, ignore.case = TRUE)] <- "Brain"
df$Tissue[grep("Dermal", df$Name, ignore.case = TRUE)] <- "Skin"
df$Tissue[grep("Ovary", df$Name, ignore.case = TRUE)] <- "Reproductive"
df$Tissue[grep("Placenta", df$Name, ignore.case = TRUE)] <- "Reproductive"
df$Tissue[grep("Stomach", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Esoph", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Rectal", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Colon", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Artery", df$Name, ignore.case = TRUE)] <- "Cardiovascular"
df$Tissue[grep("Heart", df$Name, ignore.case = TRUE)] <- "Cardiovascular"
df$Tissue[grep("hematopoi", df$Name, ignore.case = TRUE)] <- "lmmune or blood"
df$Tissue[grep("Ventricle", df$Name, ignore.case = TRUE)] <- "Cardiovascular"
df$Tissue[grep("Gastric", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Vagina", df$Name, ignore.case = TRUE)] <- "Reproductive"
df$Tissue[grep("Testis", df$Name, ignore.case = TRUE)] <- "Reproductive"
df$Tissue[grep("Intestine", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Aorta", df$Name, ignore.case = TRUE)] <- "Cardiovascular"
df$Tissue[grep("Atrium", df$Name, ignore.case = TRUE)] <- "Cardiovascular"
df$Tissue[grep("Nerve", df$Name, ignore.case = TRUE)] <- "Brain"
df$Tissue[grep("T_helper", df$Name, ignore.case = TRUE)] <- "lmmune or blood"
df$Tissue[grep("Uterus", df$Name, ignore.case = TRUE)] <- "Reproductive"
df$Tissue[grep("Ileum", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Duodenum", df$Name, ignore.case = TRUE)] <- "Gastrointestinal"
df$Tissue[grep("Nerve", df$Name, ignore.case = TRUE)] <- "Brain"
df$Tissue[grep("Nerve", df$Name, ignore.case = TRUE)] <- "Brain"
df$Tissue[grep("Pancrea", df$Name, ignore.case = TRUE)] <- "Endocrine"
df$Tissue[grep("Adrenal", df$Name, ignore.case = TRUE)] <- "Endocrine"
df$Tissue[grep("Thyroid", df$Name, ignore.case = TRUE)] <- "Endocrine"

# Write df for supplementary data:
write_tsv(df, "mtag_dyslexia_Multi_tissue_chromatin_merged.cell_type_results.txt")


#---------------------------------------------------------------------------------------------------#
# 5.2 Plot - Trait 1 Dyslexia
#---------------------------------------------------------------------------------------------------#
# Subset brain related ones only
df2 <- filter(df, Tissue == "Brain")

df2$Group <- "Adult brain"
df2$Group[grep("fetal", df2$Name, ignore.case = TRUE)] <- "Fetal brain" 
df2$Group[grep("primary_cultured_neurospheres", df2$Name, ignore.case = TRUE)] <- "Cultured"

significant_chromatin = 0.05/length(unique(df$Name))

# make Group a factor so we can use a non-alphabetic sort order
df2$Group <- factor(df2$Group, levels=c("Fetal brain", "Adult brain", "Cultured"))

# sort by Group and p-value, then add a new column with the sort orders
df2_sorted <- df2 %>% arrange(Group, Coefficient_P_value) %>% rowid_to_column(var = "sortorder")


# Plot
# use the `sortorder` column to sort the plot
chromatinplot <- ggplot(df2_sorted, aes(x=reorder(Name, sortorder), y=-log(Coefficient_P_value, 10), fill=Group)) +
  geom_bar(stat="identity", width=.75) +
  scale_y_continuous(limits = c(0,max(-log(df2$Coefficient_P_value, 10))), expand = c(0, 0)) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  geom_hline(yintercept=-log(significant_chromatin,10), color="grey60", linetype = "longdash") +
  theme_classic() + 
  theme(strip.background = element_blank()) +
  theme(text = element_text(size=10), axis.title = element_text(size = 12), strip.text = element_text(size = 10), axis.text.y = element_text(size = 8, vjust = 0), 
        axis.text.x = element_text(size = 10, angle = 90,  hjust = 1, vjust = 0.5)) +
  labs(x="", y="-log10(Pvalue)")

chromatinplot

ggsave("mtag_dyslexia_partitionedh2.png", chromatinplot, bg="white", height = 5, width = 8, scale = 2.2)
   
#---------------------------------------------------------------------------------------------------#
# DONE
#---------------------------------------------------------------------------------------------------#
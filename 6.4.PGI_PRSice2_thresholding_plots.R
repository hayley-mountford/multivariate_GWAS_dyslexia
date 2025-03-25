# Script for SNP P value threshold plots for PGIs from PRSice2
library(tidyverse)
library(viridis)

# setwd

# Based on code from https://choishingwan.github.io/PRS-Tutorial/plink_visual/

#---------------------------------------------------------------------------------------------------#
# 6.2 Age 7
#---------------------------------------------------------------------------------------------------#
# Age 7 PRS - best is 0.005
prs.result <- read_table ("prs.dyslexia.score.age7.prsice")

# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                     prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) &
                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)


# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
# PLUS VIRIDIS FOR BETTER COLOUR SCALE
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
#  Specify that we want to print p-value on top of the bars
   geom_text(
    aes(label = paste(print.p)),
    vjust = -0.5,
    hjust = 0,
    angle = 60,
    cex = 4,
    parse = T)  +
  # Add * for most thresholded value at 0.005
  annotate("text", x = 11.9, y = 0.025, label = "*", size = unit(15, "pt")) +
  # Specify the range of the plot, *1.25 to provide enough space for the p-values
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  # Specify the axis labels
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PGI model fit:  ", R ^ 2))) +
  # Draw a bar plot
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  # Specify the colors
  scale_fill_viridis()+
  # Some beautification of the plot
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save the plot
ggsave("dyslexia_PRS_in_NCDS_age7.bar.png", height = 7, width = 10)

#---------------------------------------------------------------------------------------------------#
# 6.2 Age 11
#---------------------------------------------------------------------------------------------------#
## Age 11 PRS - best is 0.1
prs.result <- read_table ("prs.dyslexia.score.age11.prsice")

# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                     prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) &
                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)

ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_text(
    aes(label = paste(print.p)),
    vjust = -0.5,
    hjust = 0,
    angle = 60,
    cex = 4,
    parse = T)  +
  # Add * for most thresholded value
  annotate("text", x = 12.9, y = 0.025, label = "*", size = unit(15, "pt")) +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PGI model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_viridis()+
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save the plot
ggsave("dyslexia_PRS_in_NCDS_age11.bar.png", height = 7, width = 10)

#---------------------------------------------------------------------------------------------------#
# 6.2 Age 16
#---------------------------------------------------------------------------------------------------#
## Age 16 PRS - best is 0.005
prs.result <- read_table ("prs.dyslexia.score.age16.prsice")

# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                     prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) &
                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)

ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_text(
    aes(label = paste(print.p)),
    vjust = -0.5,
    hjust = 0,
    angle = 60,
    cex = 4,
    parse = T)  +
  # Add * for most thresholded value
  annotate("text", x = 11.9, y = 0.018, label = "*", size = unit(15, "pt")) +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PGI model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_viridis()+
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save the plot
ggsave("dyslexia_PRS_in_NCDS_age16.bar.png", height = 7, width = 10)


#---------------------------------------------------------------------------------------------------#
# 6.2 Age 23
#---------------------------------------------------------------------------------------------------#
## Age 23 PRS - best is 1
prs.result <- read_table ("prs.dyslexia.score.age23.prsice")

# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                     prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) &
                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)

ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_text(
    aes(label = paste(print.p)),
    vjust = -0.5,
    hjust = 0,
    angle = 60,
    cex = 4,
    parse = T)  +
  # Add * for most thresholded value
  annotate("text", x = 16.9, y = 0.04, label = "*", size = unit(15, "pt")) +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PGI model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_viridis()+
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save the plot
ggsave("dyslexia_PRS_in_NCDS_age23.bar.png", height = 7, width = 10)


#---------------------------------------------------------------------------------------------------#
# 6.2 Age 33
#---------------------------------------------------------------------------------------------------#
## Age 33 PRS - best is 0.001
prs.result <- read_table ("prs.dyslexia.score.age33.prsice")

# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                     prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) &
                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)

ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_text(
    aes(label = paste(print.p)),
    vjust = -0.5,
    hjust = 0,
    angle = 60,
    cex = 4,
    parse = T)  +
  # Add * for most thresholded value
  annotate("text", x = 10.9, y = 0.026, label = "*", size = unit(15, "pt")) +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PGI model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_viridis()+
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save the plot
ggsave("dyslexia_PRS_in_NCDS_age33.bar.png", height = 7, width = 10)

#---------------------------------------------------------------------------------------------------#
# 6.2 All ages reading composite
#---------------------------------------------------------------------------------------------------#
## Reading composite PRS - best is 0.005
prs.result <- read_table ("prs.dyslexia.score.all.prsice")

# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                     prs.result$print.p == 0] <-
  format(prs.result$P[!is.na(prs.result$print.p) &
                        prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)

ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  geom_text(
    aes(label = paste(print.p)),
    vjust = -0.5,
    hjust = 0,
    angle = 60,
    cex = 4,
    parse = T)  +
  # Add * for most thresholded value
  annotate("text", x = 11.9, y = 0.027, label = "*", size = unit(15, "pt")) +
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PGI model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  scale_fill_viridis()+
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save the plot
ggsave("dyslexia_PRS_in_NCDS_composite.bar.png", height = 7, width = 10)

#---------------------------------------------------------------------------------------------------#
# DONE
#---------------------------------------------------------------------------------------------------#

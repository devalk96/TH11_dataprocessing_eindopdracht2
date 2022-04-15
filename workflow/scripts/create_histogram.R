#!/usr/bin/env Rscript

## ---------------------------
## Creates a custom histogram for interpretation of predicted heart disease rate with age
##
## Author: Sander Bouwman
## Date Created: 02-04-2022
## Copyright (c) Sander Bouwman, 2022
## ---------------------------

library(ggplot2)

# Get commandline args
args = commandArgs(trailingOnly=TRUE)

input_path <- args[1]
output_path <- args[2]

# Prevent pdf file from being rendered as it is not needed
pdf(NULL)

# Render histogram plot
heartdata <- read.csv(input_path)
ggplot(heartdata, aes(Age, fill=HeartDisease)) +
  geom_histogram(bins=12) +
  geom_vline(aes(xintercept=mean(Age[HeartDisease=="False"]),
             colour="Mean age with HD"),
             size=1, linetype="dashed") +
  scale_color_manual(name="Means", values =c("dodgerblue1", "red")) +
    guides(col = guide_legend(override.aes = list(shape = 15, size=5))) +
  ggtitle("Histogram of age and HD")

# Save plot
ggsave(output_path)
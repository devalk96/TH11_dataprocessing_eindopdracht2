#!/usr/bin/env Rscript

## ---------------------------
## Cleans ARFF file removing question marks as these mean missing values.
##
## Author: Sander Bouwman
## Date Created: 14-04-2022
## Copyright (c) Sander Bouwman, 2022
## ---------------------------


library(RWeka)

args = commandArgs(trailingOnly=TRUE)

input_path <- args[1]
output_path <- args[2]

# Remove all rows containing question marks. As these are missing values
df <- read.arff(input_path)
is.na(df) <- df == "?"
df.cleaned <- na.omit(df)
write.arff(df.cleaned, output_path)

library(edgeR)
library(tidyverse)

rawCounts <- read_tsv('./analysis/4_atacseq/normalization/rawCounts.txt')

rawCounts.matrix <- as.matrix(rawCounts[,4:6])

## edgeR:: calcNormFactors
NormFactor <- calcNormFactors(object = rawCounts.matrix, method = "TMM")

## if you prefer to use the DESeq2 strategy use method="RLE" instead

## raw library size:
LibSize <- colSums(rawCounts.matrix)

## calculate size factors:
SizeFactors <- NormFactor * LibSize / 1000000

## Reciprocal, please read section below:   
SizeFactors.Reciprocal <- 1/SizeFactors

write.table(as.data.frame(SizeFactors.Reciprocal),'./analysis/4_atacseq/normalization/ATACseq_ScaleFactors_08.05.21.txt',sep = '\t')

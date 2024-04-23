# DESeq2 experimental design and interpretation
# https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html


# load packages
library(tidyverse)
library(DESeq2)
library(knitr)


# 1. Example 1: two-group comparison
?makeExampleDESeqDataSet

# First make some example data.
dds <- makeExampleDESeqDataSet(n=10000,m=6)
assay(dds)[ 1:10,]

# This is a very simple experiment design with two conditions.
colData(dds)

dds <- DESeq(dds)
resultsNames(dds)

# explore results, here 'A' is the control or reference group
res <- results(dds, contrast=c("condition","B","A"))
# arrange the padj column from smallest to largest values/ascending order.
res <- res[order(res$padj),]

# Subset the res data frame with first 5 rows and exclude columns 3 and 4, then format using kable
library(knitr)
kable(res[1:5,-(3:4)])

# use 'B' as control or reference group
res <- results(dds, contrast=c("condition","A","B"))
# Find the index of the minimum padj value
ix = which.min(res$padj)
res <- res[order(res$padj),]
kable(res[1:5,-(3:4)])

# most significant gene.
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )






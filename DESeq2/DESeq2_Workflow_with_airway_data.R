# Thai script from Jubayer Hossain

# Install required
BiocManager::install("airway")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")

# load packages
library(tidyverse)
library(airway)
library(DESeq2)
library(EnhancedVolcano)


# get info about data
?airway

# get data
# 1. count table ~ expression level
data(airway)
airway
class(airway)

# extract count data
counts_data <- assay(airway)
head(counts_data)

# 2. metadata ~ sample information
col_data <- as.data.frame(colData(airway))

# Step 1. Preparing data
col_data <- col_data |>
  select(c(2,3)) |>
  rename(dexamethasone = dex) |>
  rename(cellline = cell)

col_data$dexamethasone <- gsub('trt', 'treated', col_data$dexamethasone)
col_data$dexamethasone <- gsub('untrt', 'untreated', col_data$dexamethasone)
col_data


# making sure the row names in col_data matches to column names in counts_data
all(colnames(counts_data) %in% rownames(col_data))

# are they in the same order?
all(colnames(counts_data) == rownames(col_data))


# Step 2: construct a DESeqDataSetFromMatrix data object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~dexamethasone)



# pre-filtering: removing rows with low gene counts
# keep rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# reference category ~ set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")


# Step 3. Run DESeq2
dds <- DESeq(dds)

# 1. Estimating Size Factors
# Purpose: Size factors account for differences in library size (sequencing depth) between samples.

# 2. Estimating Dispersions
# Purpose: Dispersion estimates capture the variability in gene expression not accounted for by the size factors.

# 3. Gene-wise Dispersion Estimates
# Purpose: Obtain gene-specific dispersion estimates to better model expression variability.

# 4. Mean-Dispersion Relationship
# Purpose: Examine the relationship between the mean expression and dispersion, which helps in understanding the data and choosing appropriate statistical models.

# 5. Final Dispersion Estimates
# Purpose: Based on the mean-dispersion relationship, obtain final dispersion estimates.

# 6. Fitting Model and Testing
# Purpose: Fit the `negative binomial model` to the data and perform hypothesis testing for `differential expression`.

# Why negative binomial model?
# 1. The negative binomial (NB) distribution is commonly used in the analysis of `count data`, such as RNA-seq data, due to its ability to model `overdispersion.`
# 2. Flexibility for Overdispersion
# 3. Biological Variability
# 4. Better Fit to Real Data
# 5. Statistical Inference

# save as result
res <- results(dds)
res

# exploring results
summary(res)

# working with alpha (significance level)
# define alpha (significance level)
# alpha = 0.01 ~ 99% CI
# alpha = 0.05 ~ 95% CI

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

resultsNames(dds)

# contrast
contrast_res <- results(dds, contrast = c("dexamethasone", "treated", "untreated"))
contrast_res


# visualizations (What? Why? How?)

# MA Plot (M vs A Plot)
# What: An MA plot represents the log-fold change (M) on the y-axis and the average expression (A) on the x-axis for each gene or feature.
# Why: It is used to visualize differential expression between two conditions.The log-fold change (M) gives an idea of the magnitude of change, and the average expression (A) helps identify if the change is dependent on the expression level.

# Interpretations
# 1. Points on the plot represent genes.
# 2. Genes with significant differential expression are often `found at the extremes` of the plot.
# 3. Upregulated genes are at the top
# 4. Downregulated genes at the bottom
# 5. Non-differentially expressed genes are centered around zero on the y-axis.

plotMA(dds)

dds <- makeExampleDESeqDataSet(n=10000,m=18)


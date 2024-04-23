# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Check Bioconductor version
BiocManager::version()

# Check available packages
BiocManager::available()

# Install single Bioconductor package
BiocManager::install("Package+_Name")

# Install multiple Bioconductor package
BiocManager::install(c("Package_1", "Package_2", "Package_3"))

# Load packages
library(package_name)

# Package version
packageVersion("GenomicRanges")

# check Reference manual
?GenomicRanges

# Check Vignettes
browseVignettes("GenomicRanges")


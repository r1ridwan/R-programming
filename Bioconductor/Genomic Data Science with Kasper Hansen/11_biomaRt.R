library(biomaRt)
BiocManager::install("biomaRt")

#Package Name: biomaRt
#Description: Interface to BioMart databases (i.e. Ensembl)
#Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/biomaRt.html
#Documentation: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

## ----listMarts-----------------------------------------------------------
head(listMarts())
mart <- useMart("ensembl")
mart
listDatasets(mart)
head(listDatasets(mart))
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
ensembl

## ----getBMex-------------------------------------------------------------
values <- c("202763_at","209310_s_at","207500_at")
getBM(attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"),
      filters = "affy_hg_u133_plus_2", values = values, mart = ensembl)

## ----listAttributes------------------------------------------------------
attributes <- listAttributes(ensembl)
head(attributes)
nrow(attributes)
tail(attributes)

#---------ListFilters------------------------------
filters <- listFilters(ensembl)
head(filters)
nrow(filters)

## ----listPages-----------------------------------------------------------
attributePages(ensembl)
attributes <- listAttributes(ensembl, page = "feature_page")
head(attributes)
nrow(attributes)
attributes <- listAttributes(ensembl, page = "sequences")
attributes
nrow(attributes)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


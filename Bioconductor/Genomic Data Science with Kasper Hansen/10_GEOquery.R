library(GEOquery)

#Package Name: GEOquery
#Description: Get data from NCBI Gene Expression Omnibus (GEO)
#Bioconductor URL: https://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
#Documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html

## ----getData-------------------------------------------------------------
#Download data from GEO with the GSE11675 Accession number
eList <- getGEO("GSE11675")
class(eList)
length(eList)
names(eList)
eData <- eList[[1]]
eData


## ----pData---------------------------------------------------------------
names(pData(eData))
eData$geo_accession
eData$submission_date
eData$type

## ----getGEOsupp----------------------------------------------------------
eList2 <- getGEOSuppFiles("GSE11675")
eList2
tarArchive <- rownames(eList2)[1]
tarArchive

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


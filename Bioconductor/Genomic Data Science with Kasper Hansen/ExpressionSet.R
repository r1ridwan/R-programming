library(Biobase)
library(ALL)
library(hgu95av2.db)

#Package Name: ALL
#Description: A data package- The data consist of microarrays from 128 different individuals with acute lymphoblastic leukemia (ALL)
#Bioconductor URL: https://bioconductor.org/packages/release/data/experiment/html/ALL.html
#Documentation: PDF


#Load data from ALL data package
data(ALL)
ALL

## ----help, eval=FALSE----------------------------------------------------
?ALL

## ----experimentData------------------------------------------------------
experimentData(ALL)

#See the abstract of this data package
abstract(ALL)

## ----exprs---------------------------------------------------------------
exprs(ALL)[1:4, 1:4]

## ----names---------------------------------------------------------------
sampleNames(ALL)
head(sampleNames(ALL))
featureNames(ALL)
head(featureNames(ALL))

## ----phenotypic-Data(pData)---------------------------------------------------------------
head(pData(ALL))

## ----dollar--------------------------------------------------------------
pData(ALL)$sex
head(pData(ALL)$sex)
ALL$sex
head(ALL$sex)


## ----subset--------------------------------------------------------------
#Subset first five sample
ALL[,1:5]
#Subset first ten features
ALL[1:10,]
#Subset first ten features and first five samples
ALL[1:10,1:5]

## ----subset2-------------------------------------------------------------
#Subset sample 3, 2, 1
ALL[, c(3,2,1)]
#Subset sex of sample 1, 2, 3
ALL$sex[c(1,2,3)]
#Subset sex of sample 3, 2, 1
ALL[, c(3,2,1)]$sex

## ----featureData---------------------------------------------------------
featureData(ALL)

## ----annotation----------------------------------------------------------
ids <- featureNames(ALL)[1:5]
ids

## ----annotation2---------------------------------------------------------
?hgu95av2.db
library(hgu95av2.db)
as.list(hgu95av2ENTREZID[ids])

## ----varLabels-----------------------------------------------------------
pD = phenoData(ALL)
pD
varLabels(pD)
names(pData(ALL))
pData(ALL)
pData(phenoData(ALL))

## ----varLabels2----------------------------------------------------------
#change column 2 name
varLabels(pD)[2] <- "Age at diagnosis"
pD
colnames(pD)[1:3]
varLabels(pD)[1:3]
#Change column 3 name
varLabels(pD)[3] <- "SEX"
varLabels(pD)[1:3]


## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


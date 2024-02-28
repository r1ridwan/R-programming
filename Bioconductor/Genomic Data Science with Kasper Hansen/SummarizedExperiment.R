library(GenomicRanges)
install.packages("airway")
library(airway)


#Package Name: airway
#Description: RangedSummarizedExperiment for RNA-Seq in airway smooth muscle cells, by Himes et al PLoS One 2014
#Bioconductor URL: https://bioconductor.org/packages/release/data/experiment/html/airway.html
#Documentation: https://bioconductor.org/packages/release/data/experiment/vignettes/airway/inst/doc/airway.html
#SummerrizedExperiment: https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html

#Load data from airway package
data(airway)
airway

## ----colData-------------------------------------------------------------
colData(airway)

## ----getColumn-----------------------------------------------------------
airway$cell
airway$SampleName
airway$dex

## ----exptData------------------------------------------------------------
exptData(airway)

## ----names---------------------------------------------------------------
colnames(airway)
head(rownames(airway))

## ----assay---------------------------------------------------------------
airway
assayNames(airway)
assays(airway)
head(assay(airway, "counts"))[1:4, 1:4]

## ----rowRanges-----------------------------------------------------------
rowRanges(airway)
length(rowRanges(airway))
dim(airway)

## ----numberOfExons-------------------------------------------------------
length(rowRanges(airway))
sum(elementLengths(rowRanges(airway)))

## ----start---------------------------------------------------------------
start(rowRanges(airway))
start(airway)

## ----subsetByOverlaps----------------------------------------------------
gr <- GRanges(seqnames = "1", ranges = IRanges(start = 1, end = 10^7))
subsetByOverlaps(airway, gr)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


library(GenomicRanges)
library(IRanges)
library(plyranges)
library(GenomeInfoDb)
library(AnnotationHub)
library(devtools)
library(rtracklayer)
library(dplyr)
library(BiocManager)
BiocManager::install("AnnotationHub")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

# Package Name: AnnotationHub
# Description: Client to access AnnotationHub resources
# Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html
# Access the AnnotationHub Web Service:https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub.html
# AnnotationHub How-To's: https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub-HOWTO.html
# Troubleshoot The Hubs: https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/TroubleshootingTheCache.html
# Kasper Hansen: http://kasperdanielhansen.github.io/genbioconductor/html/AnnotationHub.html


## ----annoHub,results="hide"----------------------------------------------
ah = AnnotationHub()
ah
length(ah) 

## ----no1-----------------------------------------------------------------
ah[1]
ah[2]

## ----look----------------------------------------------------------------
unique(ah$dataprovider)
unique(ah$rdataclass)
unique(ah$species)
head(unique(ah$species))

## ----subset--------------------------------------------------------------
ah <- subset(ah, species == "danio rerio")
ah

## ----query---------------------------------------------------------------
query(ah, "H3K4me3")
query(ah, c("H3K4me3", "Gm12878"))
dm <- query(ah, c("ChainFile", "UCSC", "Drosophila melanogaster"))
dm
ahs <- query(ah, c("inparanoid8", "ailuropoda"))
ahs

#--------Retrieve-metadata-------------------------------------------
df <- mcols(dm)
df
length(ah)
dm ["AH15146"]
dm[["AH15146"]]


#------------Display-data------------------------------------------
# It will open shiny we application and then select your data rows and click send
# Then execute the function, it will open the selected data rows
ah2 <- display(ah)
ah2

#---------------configuration------------------------------------
ah
snapshotDate(ah)
pd <-possibleDates(ah)
pd
snapshotDate(ah) <- pd[1]
snapshotDate(ah)
ah



#------------------------------------------------------------------------------------------
#URL: http://kasperdanielhansen.github.io/genbioconductor/html/Usecase_AnnotationHub_GRanges.html

## ----ahub_species--------------------------------------------------------
ah <- AnnotationHub()
ah
ah <- subset(ah, species == "Homo sapiens")
ah


## ----ahub_histone-protein-------------------------------------------------------
qhs <- query(ah, "H3K4me3")
qhs <- query(qhs, "Gm12878")

## ----ahub_look-----------------------------------------------------------
qhs
qhs[2]
qhs[4]
qhs[6]

## ----ahub_closerlook-----------------------------------------------------
qhs$title
qhs$dataprovider

## ----ahub_twoGR----------------------------------------------------------
gr1 <- qhs[[1]]
gr1
gr2 <- qhs[[2]]
gr2
gr3 <- subset(qhs, title == "wgEncodeUwHistoneGm12878H3k4me3StdPkRep1.narrowPeak.gz")[[1]]
gr3
gr2 <- subset(qhs, title == "E116-H3K4me3.narrowPeak.gz")[[1]]
gr2

## ----ahub_summary--------------------------------------------------------
summary(width(gr1))
summary(width(gr2))
table(width(gr1))
table(width(gr2))

## ----ahub_refseq---------------------------------------------------------
qhs <- query(ah, "RefSeq")
qhs

## ----ahub_refseq_genome--------------------------------------------------
qhs$genome

## ----ahub_histone_genome-------------------------------------------------
genome(gr1)
genome(gr2)

#-------------get-genes-info-----------------------------------------------
genes <- qhs[[1]]
genes
peaks = gr1
peaks

## ----ahub_get_refseq-----------------------------------------------------
refseq <- qhs[qhs$genome == "hg19" & qhs$title == "RefSeq Genes"]
refseq
refseq <- refseq[[1]] ## Downloads

## ----refseq--------------------------------------------------------------
refseq

## ----ahub_refseq_name----------------------------------------------------
table(table(refseq$name))

## ----promoters-----------------------------------------------------------
prom <- promoters(refseq)
prom
table(width(prom))
args(prom)


## ----findOverlaps--------------------------------------------------------
ov <- findOverlaps(promoters, peaks)
ov

## ----queryHits-----------------------------------------------------------
length(peaks)
length(prom)
length(unique(queryHits(ov)))
length(unique(subjectHits(ov)))
length(unique(queryHits(ov))) / length(peaks)
length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE))
length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE)) / length(peaks)
length(subsetByOverlaps(prom, peaks, ignore.strand = TRUE)) / length(prom)

## ----subjectHits---------------------------------------------------------
length(unique(subjectHits(ov)))
length(unique(subjectHits(ov))) / length(promoters)


## ----widthPercentage-----------------------------------------------------
sum(width(reduce(gr1))) / 10^6
sum(width(reduce(prom))) / 10^6

## ----size----------------------------------------------------------------
sum(width(intersect(gr1, prom))) / 10^6

## ----size2---------------------------------------------------------------
sum(width(intersect(gr1, prom, ignore.strand = TRUE))) / 10^6

## ----widthPercentage2----------------------------------------------------
sum(width(reduce(prom))) / 10^6
sum(width(reduce(prom, ignore.strand = TRUE))) / 10^6

## ----promInOut-----------------------------------------------------------
prom <- reduce(prom, ignore.strand = TRUE)
peaks <- reduce(gr1)
both <- intersect(prom, peaks)
only.prom <- setdiff(prom, both)
only.peaks <- setdiff(peaks, both)
overlapMat <- matrix(0,, ncol = 2, nrow = 2)
colnames(overlapMat) <- c("in.peaks", "out.peaks")
rownames(overlapMat) <- c("in.promoters", "out.promoter")
overlapMat[1,1] <- sum(width(both))
overlapMat[1,2] <- sum(width(only.prom))
overlapMat[2,1] <- sum(width(only.peaks))
overlapMat[2,2] <- 3*10^9 - sum(overlapMat)
round(overlapMat / 10^6, 2)

## ----oddsratio-----------------------------------------------------------
oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
oddsRatio

## ----oddsRatio2----------------------------------------------------------
overlapMat[2,2] <- 1.5*10^9
oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
oddsRatio

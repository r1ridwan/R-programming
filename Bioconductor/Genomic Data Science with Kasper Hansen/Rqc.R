library(ShortRead)
library(Rqc)


#Package Name: Rqc
#Description: Quality Control Tool for High-Throughput Sequencing Data
#Bioconductor URL: https://www.bioconductor.org/packages/release/bioc/html/Rqc.html
#Documentation URL: https://www.bioconductor.org/packages/release/bioc/vignettes/Rqc/inst/doc/Rqc.html

##============================Rqc-====================================================
#===========================================================
## Multiple and Parallal Sequence Quality
BiocManager::install("Rqc")
folder <- system.file(package="ShortRead", "fasta/SRR1971253_pass.fastq")

fastqDir <- system.file(package="ShortRead", "extdata/E-MTAB-1147")
files <- list.files(fastqDir, "fastq.gz", full.names=TRUE)
qarqc <- rqcQA(files, workers = 1)
class(qarqc)
names(qarqc)
reportFile <- rqcReport(qarqc)
browseURL(reportFile)
#-------------
filess <- "fasta/fastqqqfile.fastq"
qa <- rqcQA(filess, workers = 1)
perFileInformation(qa)
#Sample of sequenc
set.seed(111)
qa_sample <- rqcQA(filess, workers = 1, sample = TRUE, n=500)
qa_sample
#first two files are one pair and second two files are the next pari
qa_paired <- rqcQA(qa_sample, workers = 4, pair = c(1, 1, 2, 2))

#quality control html report
reportfile <- rqcReport(qa, templateFile = "myReport.Rmd")
browseURL(reportfile)
#The class of qa is rqcResultSet
methods(class = "RqcResultSet")

## Some rqc plot functions
#rqcCycleAverageQualityPcaPlot()
#rqcCycleAverageQualityPlot()
#rqcCycleBaseCallsLinePlot()
#rqcCycleBaseCallsPlot()
#rqcCycleGCPlot()
#rqcCycleQualityBoxPlot()
#rqcGroupCycleAverageQualityPlot()
#rqcReadFrequencyPlot()
#rqcReadQualityBoxPlot()
#rqcReadQualityPlot()
#rqcReadWidthPlot()
#rqcCycleQualityPlot()

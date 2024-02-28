library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

#Package Name: rtracklayer
#Description: R interface to genome annotation files and the UCSC genome browser
#Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/rtracklayer.html
#Documentation: PDF

?import
?BigWigFile

## ----AnnotationHub----------------------------------------------------------------
ahub <- AnnotationHub()
table(ahub$rdataclass)

## ----granges-------------------------------------------------------------
#Load GRanges rdataclass data from annotationhub
ahub.gr <- subset(ahub, rdataclass == "GRanges" & species == "Homo sapiens")
gr <- ahub.gr[[1]]
gr
seqinfo(gr)

## ----BigWig--------------------------------------------------------------
#Load BigWig rdataclass file data from annotationhub
ahub.bw <- subset(ahub, rdataclass == "BigWigFile" & species == "Homo sapiens")
ahub.bw
bw <- ahub.bw[[1]]
bw

#Import BigWig file
gr.chr22 <- import(bw, which=GRanges("chr22", ranges= IRanges(1, 10^8)))
gr.chr22
rle.chr22 <- import(bw, which=GRanges("chr22", ranges= IRanges(1, 10^8)), as="Rle")
rle.chr22
rle.chr22$chr22


## ----importBigWig3-------------------------------------------------------
gr.chr22 <- GRanges(seqnames = "chr22",
                    ranges = IRanges(start = 1, end = seqlengths(gr)["chr22"]))
out.chr22 <- import(bw, which = gr.chr22, as = "Rle")
out.chr22[["chr22"]]

## ----liftOver------------------------------------------------------------
ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
query(ahub.chain, c("hg18", "hg19"))
chain <- ahub.chain[ahub.chain$title == "hg19ToHg18.over.chain.gz"]
chain <- chain[[1]]
gr.hg18 <- liftOver(gr, chain)
gr.hg18

## ----liftOver2-----------------------------------------------------------
table(elementLengths(gr.hg18))

## ----tabixIndex----------------------------------------------------------
library(Rsamtools)
from <- system.file("extdata", "ex1.sam", package="Rsamtools",
                    mustWork=TRUE)
from
to <- tempfile()
zipped <- bgzip(from, to)
idx <- indexTabix(zipped, "sam")


#Load Microarray data from rtracklayer--------------------
data("targets")
head(targets)
targetranges <- IRanges(targets$start, targets$end)
targetstrack <- with(targets, GRangesForUCSCGenome("hg18", chrom, targetranges, strand, name, target))
genome(targetstrack)
head(seqlengths(targetstrack))

#Accessing Track informtions--------------------------------
head(seqnames(targetstrack))
head(start(targetstrack))

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


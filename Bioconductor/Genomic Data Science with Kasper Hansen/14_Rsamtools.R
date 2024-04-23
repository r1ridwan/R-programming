library(Rsamtools)

#Package Name: Rsamtools
#Description: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import
#Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/Rsamtools.html
#Documentation URL: PDF


## ----bamPath-------------------------------------------------------------
bamPath <- system.file("extdata", "ex1.bam", package="Rsamtools")
bamPath
bamFile <- BamFile(bamPath)
bamFile

## ----bamFileInfo---------------------------------------------------------
seqinfo(bamFile)

## ----scanBam-------------------------------------------------------------
aln <- scanBam(bamFile)
length(aln)
class(aln)

## ----lookAtBam-----------------------------------------------------------
aln <- aln[[1]]
names(aln)
lapply(aln, function(xx) xx[1])

## ----yieldSize-----------------------------------------------------------
yieldSize(bamFile) <- 1
open(bamFile)
scanBam(bamFile)[[1]]$seq
scanBam(bamFile)[[1]]$seq
## Cleanup
close(bamFile)
yieldSize(bamFile) <- NA

## ----ScanBamParams-------------------------------------------------------
gr <- GRanges(seqnames = "seq2",
              ranges = IRanges(start = c(100, 1000), end = c(1500,2000)))
params <- ScanBamParam(which = gr, what = scanBamWhat())
aln <- scanBam(bamFile, param = params)
names(aln)
head(aln[[1]]$pos)

## ----summary-------------------------------------------------------------
quickBamFlagSummary(bamFile)

## ----BamViews------------------------------------------------------------
bamView <- BamViews(bamPath)
aln <- scanBam(bamView)
names(aln)

## ----BamViews2-----------------------------------------------------------
bamRanges(bamView) <- gr
aln <- scanBam(bamView)
names(aln)
names(aln[[1]])

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


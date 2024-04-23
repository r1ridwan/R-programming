library(GenomicRanges)
library(IRanges)
library(plyranges)
library(GenomeInfoDb)
library(AnnotationHub)
library(devtools)
library(rtracklayer)
library(BSgenome)

#Package Name: GenomicRanges: data structure for storing genomic intervals
#Description: Representation and manipulation of genomic intervals
#Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
#Vignettes PDF
#Documentation URL: https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
#Kasper Hansen: http://kasperdanielhansen.github.io/genbioconductor/


# Create-GRanges-Object---------------
gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),
              ranges = IRanges(start = c(1,3,5), width = 3))
     
gr 

gr2 <- GRanges(seqnames = "chr1", 
               strand = c("+", "-", "*"),
               ranges = IRanges(start = c(1, 3, 5), width = 3),
               gene_id = c(10025, 10026, 10027),
               score = c(10, 15, 12))
gr2

# (+) Indicate Forward strand
# (-) Indicate Reverse strand
# (*) Indicate either unknown strand or indicate presence in both strand

gr3 <- GRanges(seqnames = c("chr1", "chr2", "chr3"), 
               strand = c("+", "-", "*"),
               ranges = IRanges(start = c(1, 3, 5), width = 3),
               gene_id = c(10025, 10026, 10027),
               score = c(10, 15, 12))
gr3


#----------GenomicRanges--------------------------------------------------------
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer2")
library(BSgenome.Scerevisiae.UCSC.sacCer2)
yeast <- BSgenome.Scerevisiae.UCSC.sacCer2
seqnames(yeast)
chrI <- yeast$chrI
(myGR <- GRanges("chrI:200-300"))
# Get GenomicRanges Accessor
methods(class = "GRanges")
# Chromosome Nmaes
seqnames(myGR)
# IRanges object for ranges
ranges(myGR)
# Metadata columns
mcols(myGR)
# Sequence Information
seqinfo(myGR)
# Get Genome Name
genome(myGR)

## ----flank---------------------------------------------------------------------------
flank(gr, 2) # Positive strand value shift left(decrease) and negative strand value shift right(increase)
flank(gr, 3, start = FALSE) # Psitive strand value shift right (increase) and negative value shift left (decrease)

## ----------seqinfo-seqlengths-------------------------------------------------------------
seqinfo(gr2)
seqlengths(gr2) <- c("chr1" = 10)
seqinfo(gr)
seqlevels(gr)
seqlengths(gr)
seqlengths(gr3) <- c(249250621, 243199373, 198022430)
seqlengths(gr3)

#----------------ranges----------------------------------------
ranges(gr3)

#-------------------strand------------------------------------
strand(gr3)

## ----gaps----------------------------------------------------------------
gaps(gr3)

#--------------------names--------------------------------------------
names(gr3)

#----------------length----------------------------------------
length(gr3)

#-----------granges-------------------------------------------------
#The genomic ranges can be extracted without corresponding metadata with granges
granges(gr3)

## ----gr2-giving-new-seqnames-----------------------------------------------------------------
seqnames(gr3)
seqlevels(gr) <- c("chr1", "chr2")
seqnames(gr) <- c("chr1", "chr2", "chr1")
gr

seqlevels(gr) <- c("chr1", "chr2", "chr3")
seqnames(gr) <- c("chr1", "chr2", "chr3")
gr

## ----sorting-seqlevels----------------------------------------------------------------
sort(gr)
seqlevels(gr) <- c("chr2", "chr1")
sort(gr)

## ----genome--------------------------------------------------------------
genome(gr) <- "hg19"
gr
seqinfo(gr)

## ----gr-error, error=TRUE------------------------------------------------
# This becomes valuable when you deal with data from different genome versions, because it allows R to throw an error when you compare two GRanges from different genomes, like
gr2 <- gr
genome(gr2) <- "hg18"
findOverlaps(gr, gr2)

##------Splitting-and-Combining-Grange-Objects----------------------------------------------
sp <- split(gr4, rep(1:2, each= 5))
sp[1]
sp[2]
c(sp[[1]], sp[[2]]) #Concatenated

#---------Subsetting-GRanges-Objects----------------------------------------------------
gr4[2:3]
gr4[8]
gr4[2:3, "GC"]
gr4[4:5, "score"]
strand(gr4)
ranges(gr4)
seqinfo(gr4)
length(gr4)
seqlengths(gr4)
names(gr4)
granges(gr4)
seqnames(gr4)
seqlevels(gr4)
#Elements can also be assigned to the GRanges object. Here is an example where the second row of a GRanges object is replaced with the first row of gr4
singles <- split(gr4, names(gr4))
grMod <- gr4
grMod[2] <- singles[[1]] #2nd row is replaced by 1st row
head(grMod, n = 4)
gr4
grMod[3] <- singles[[1]] #3rd row is replaced by 1st row
head(grMod, n = 4)

grMod[1] <- singles[[10]] #1st row is replaced by 10th row
head(grMod, n = 4)

#Repeat a Row
rep(singles[[2]], times = 3) #2nd row is repeated 3 times
rep(singles[[10]], times=10) #10th ros is repeated 10 times
rep(singles[[3]], times = 12)

# Reverse the whole GRanges Objects
rev(gr4)

head(gr4, n=2)
tail(gr4, n=2)
window(gr4, start = 2, end = 4)
gr4[IRanges(start = c(2, 7), end = c(3, 9))]



#----------------------------------------------------------------------------------------------
#URL: http://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html
# ----DataFrame-----------------------------------------------------------
ir <- IRanges(start = 1:2, width = 3)
ir
df1 <- DataFrame(iranges = ir)
df1
df1$iranges
df2 <- data.frame(iranges = ir)
df2

## ----GRanges-------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),
              ranges = IRanges(start = c(1,3,5), width = 3))
values(gr) <- DataFrame(score = c(0.1, 0.5, 0.3))
values(gr)
gr

#-------------metadata-extraction------------------------------------------
#Annotations or metadata of a grange can be extracted as a DataFrame object using the mcols accessor.
mcols(gr2)
mcols(gr3)$score
mcols(gr3)$GC

## ----grdollar------------------------------------------------------------
gr$score
gr
gr$score2 = gr$score * 0.2
gr
gr$score3 = gr$score2 / 3
gr

## ----findOverlaps_setup--------------------------------------------------
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = "*",
               ranges = IRanges(start = c(1, 3, 5), width = 3))
gr2
gr

## ----findOverlaps--------------------------------------------------------
findOverlaps(gr, gr2)
findOverlaps(gr, gr2, ignore.strand = TRUE)

## ----subsetByOverlaps----------------------------------------------------
gr
subsetByOverlaps(gr, gr2)
gr2
subsetByOverlaps(gr2, gr)

## ----makeGRangesFromDataFrame--------------------------------------------
#Convert datagrame into GRanges
df <- data.frame(
  chr = "chr1", 
  start = 1:3, 
  end = 4:6, 
  score = 7:9,
  strand = c("+", "-", "+"),
  GC = c(0.25, 0.24, 0.23))
df
makeGRangesFromDataFrame(df)
makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
gr5 <- as(df, "GRanges") #transform df into GRanges
gr5

#-----------------make-granges-from-dataframe---------------------------------
#Transform dataframe into GRanges
gr <- 
  data.frame(seqnames = sample(c("chr1", "chr2", "chr3"), 7, replace = TRUE), 
             strand = sample(c("+", "-", "*"), 7, replace = TRUE),
             start = 1:7, 
             width = 10, 
             score = runif(7)) %>%
  as_granges()
gr

##----------Summarizing by group-(plyranges-package)--------------------------------------------
gr %>%
  group_by(strand) %>% 
  summarize(mean_score = mean(score))

##-----------Filtering by group------------------------------------------
gr %>%
  group_by(strand) %>% 
  filter(score == max(score))

##=============GRangesList===============================================================
##-----------------GRangesList-------------------------------------------
#Making GRangesList or Create GRangesList
#GRangesList Accessors
?GRangesList
methods(class = "GRangesList")

gr1 <- GRanges(seqnames = "chr2", ranges = IRanges(103, 106), strand = "+", score = 5L, GC = 0.45)
gr2 <- GRanges(seqnames = c("chr1", "chr1"), ranges = IRanges(c(107, 113), width = 3),
               strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
GRangesList(gr1, gr2) #GRangesList with gr1 and gr2 GRanges
grl <- GRangesList("txA" = gr1, "txB" = gr2) #GRangesList with gr1 and gr2 GRanges
grl
#Basic GRangesList accessors
start(grl)
seqnames(grl)
ranges(grl)
strand(grl)
length(grl)
names(grl)
seqlengths(grl)
elementNROWS(grl)
shift(grl, 10)

##findOverlaps
findOverlaps(grl, gr1)

#"isEmpty" tests if a GRangesList object contains anything.
isEmpty(grl)
mcols(grl) <- c("Transcript A","Transcript B")
mcols(grl)
#Element-level metadata can be retrieved by unlisting the GRangesList, and extracting the metadata
mcols(unlist(grl))

#Combining GRangesList objects--------------------------------
#GRangesList objects can be unlisted to combine the separate GRanges objects that they contain as an expanded GRanges.
ul <- unlist(grl)
ul
#Combining two grangeslist into single list
grl1 <- GRangesList(
  gr1 = GRanges("chr2", IRanges(3, 6)),
  gr2 = GRanges("chr1", IRanges(c(7,13), width = 3)))
grl2 <- GRangesList(
  gr1 = GRanges("chr2", IRanges(9, 12)),
  gr2 = GRanges("chr1", IRanges(c(25,38), width = 3)))

pc(grl1, grl2)
grl3 <- c(grl1, grl2)
regroup(grl3, names(grl3))

#------Basic interval operations for GRangesList objects------------------------
start(grl)
end(grl)
width(grl)
# sum of widths of each grl element
sum(width(grl))
shift(grl, 20)
coverage(grl)

#-------Subsetting GRangesList objects--------------------------
grl
grl[1]
grl[2]
grl[[1]]
grl["txA"]
grl$txB
grl[2, "score"]
grl["txB", "GC"]
rep(grl[[2]], times = 5)
rev(grl)
head(grl, n=1)
head(grl, n=2)
tail(grl, n=1)
window(grl, start = 1, end=1)
window(grl, start = 1, end=2)
grl[IRanges(start = 1, end = 1)]
grl[IRanges(start = 1, end = 2)]
grl[IRanges(start = 2, end = 2)]

#Looping over GRangesList objects
lapply(grl, length)
sapply(grl, length)

grl2 <- shift(grl, 10)
names(grl2) <- c("shiftTxA", "shiftTxB")
mapply(c, grl, grl2)

Map(c, grl, grl2)

endoapply(grl, rev)
mendoapply(c, grl, grl2)

Reduce(c, grl)

gr <-unlist(grl)
gr$log_score <- log(gr$score)
grl <- relist(gr, grl)
grl


#------------------------------------------------------------------------------------------
#URL: http://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_seqinfo.html
# Package URL: https://bioconductor.org/packages/devel/bioc/html/GenomeInfoDb.html
# Documentation: https://bioconductor.org/packages/devel/bioc/vignettes/GenomeInfoDb/inst/doc/GenomeInfoDb.pdf
##-------GenomeInfoDb------------------------------------------------
#--------------genomeStyles-----------------------------------------
seqmap <- genomeStyles()
head(seqmap, n = 3)
#Oragnism's supported by GenomeInfoDb can be found by :
names(genomeStyles())
#List out first five entries for Homo_sampiens
head(genomeStyles("Homo_sapiens"), 5)
#Checking a given style whether it is supported by GenomeInfoDb for a given species
"UCSC" %in% names(genomeStyles("Homo_sapiens"))

## ----seqlevelsForce-----------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
gr
seqlevels(gr4, force = TRUE) = "chr1"
gr

## ----dropSeqlevels-------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
gr
dropSeqlevels(gr4, "chr1")
gr
keepSeqlevels(gr, "chr2")
?dropSeqlevels

## ----keepStandard--------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chrU345"),
              ranges = IRanges(start = 1:2, end = 4:5))
gr
keepStandardChromosomes(gr)

## ----GRanges-------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:2, width = 2))

## ----seqStyle------------------------------------------------------------
newStyle <- mapSeqlevels(seqlevels(gr), "NCBI")
gr <- renameSeqlevels(gr, newStyle)
gr

##=================Rle=======================================================
#----------------------------GRanges-Rle---------------------------------------
#Rle stands for Run length encoding
#Rle basic Uses
rl <- Rle(c(1,1,1,1,2,2,3,3,2,2))
rl #It shows lengths and values. Here Value 1 repeated 4times, 2 repeated 2times, 3 repeated 2times, and 2 repeated 2times
runLength(rl)
runValue(rl)
as.numeric(rl)
gr4 <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1: 10,
  GC = seq(1, 0, length = 10))
gr4

##aggregate-------
ir <- IRanges(start = c(2,6), width = 2)
ir
aggregate(rl, ir, FUN = mean)
ir <- IRanges(start = c(2, 8), width = 4)
ir
vec <- as.numeric(rl)
mean(vec[2:5])
mean(vec[8:11])

##coverage-------
ir <- IRanges(start = 1:10, width = 3)
rl <- coverage(ir)
rl
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width = 3))
rl <- coverage(gr)
rl

#Views-------
vi <- Views(rl, start = c(3,7), width = 3)
vi
mean(vi)

grView <- GRanges("chr1", ranges = IRanges(start = 2, end = 7))
vi <- Views(rl, grView)
vi

vi <- Views(rl, as(grView, "RangesList"))
vi
vi[[1]]


## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()
























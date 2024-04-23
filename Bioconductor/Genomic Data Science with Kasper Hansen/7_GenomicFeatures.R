#Package Name: GenomicFeatures
#Description: Conveniently import and query gene models
#Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
#Vigenettes URL: PDF 
#Kasper Hansen: http://kasperdanielhansen.github.io/genbioconductor/

BiocManager::install("GenomicFeatures")
BiocManager::install("mirbase.db")
BiocManager::install("FDb.UCSC.tRNAs")
library(GenomicFeatures)
library(mibase.db)

#I can load pre-existing data using loadDb function 
samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                          package="GenomicFeatures")

txdb <- loadDb(samplefile)
txdb


##---------GenomicFeatures--------------------------------------------------------
#GenomicFeatures uses TXDb object to store metadata
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# Load human reference genome hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
?TxDb.Hsapiens.UCSC.hg38.knownGene
# Assign hg38 to hg, then print it
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg

# Extract all the genes in chromosome X as hg_chrXg, then print it
hg_chrXG <- genes(hg, filter = list(tx_chrom = c("chrX")))
hg_chrXG
# Extract all positive stranded genes in chromosome X, assign to hg_chrXgp, then sort it
hg_chrXgp <- genes(hg, filter = list(tx_chrom = "chrX", tx_strand = "+"))
sort(hg_chrXgp)
#If you would like to test other filters, valid names for this list are: 
#"gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand", "exon_id", "exon_name", "exon_chrom", "exon_strand", "cds_id", "cds_name", "cds_chrom", "cds_strand", and "exon_rank".
#Extract any particular gene with "gene_id"
(hg_chrXgp <- genes(hg, filter = list(tx_chrom = "chrX", gene_id = "9823")))


## ----txdb----------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
#Extract Basic quantities
genes(txdb)
transcripts(hg)
cds(hg)
exons(hg)
microRNAs(txdb)
tRNAs(txdb)
promoters(txdb)
seqlevels(txdb)
head(seqlevels(txdb))

#If you want to set chrI to be only acive, do it with the follwoing way-
seqlevels(txdb) <- "chr1"

#If you need to reset back to the original seqlevels 
seqlevels(txdb) <- seqlevels0(txdb)

#Extract quantities and group
transcriptsBy(by = c("gene", "exon", "cds"), txdb)
cdsBy(by = c("tx", "gene"), txdb)
exonsBy(by = c("tx", "gene"), txdb)
intronsByTranscript(txdb)
fiveUTRsByTranscript(txdb)
threeUTRsByTranscript(txdb)

## ----gr------------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", strand = "+", ranges = IRanges(start = 11874, end = 14409))
#Find which gene is overlapped in this granges
subsetByOverlaps(genes(txdb), gr)
subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)

#this granges has three transcripts and how they are arranged with this 6 exons. 
## ----transcripts---------------------------------------------------------
#Find which transcript is overlapped in this granges
subsetByOverlaps(transcripts(txdb), gr)

## ----exons---------------------------------------------------------------
#Find which exons are overlapped in this granges, and exons forming transcript
subsetByOverlaps(exons(txdb), gr)

## ----exonsBy-------------------------------------------------------------
#Find out which exons are forming these three transcripts on this particular granges
#Here we now finally see the structure of the three transcripts in the form of a GRangesList
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)


## ----cds-----------------------------------------------------------------
subsetByOverlaps(cds(txdb), gr)
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)["2"]

## ----transcriptLengths---------------------------------------------------
#Find out which transcript genuinely containing the coding sequence
subset(transcriptLengths(txdb, with.cds_len = TRUE), gene_id == "100287102")
sum(width(subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)[["2"]]))

##Retreiving Data using the select method------------------------------------
#how to find the UCSC transcript names that match with a set of gene IDs
seqlevels(txdb) <- "chr1"
keys <- c("100129405", "440699", "5208")
columns(txdb)
keytypes(txdb)
select(txdb, keys = keys, columns="TXNAME", keytype="GENEID")

#For the genes in the example above, find the chromosome and strand information that will go with each of the transcript names.
cols <- c("TXNAME", "TXSTRAND", "TXCHROM")
select(txdb, keys = keys, columns = cols, keytype = "GENEID")


##Methods for returning GRAnges object--------------------------------
Gr<-transcripts(txdb)
Gr[1:3]
tx_strand <- strand(Gr)
tx_strand
sum(runLength(tx_strand))
length(Gr)

#Retreive a subset of transcripts availble such as those on the positive strand of the chi1
GR <- transcripts(txdb, filter=list(tx_chrom = "chr15", tx_strand = "+"))
length(GR)
unique(strand(GR))

PR <- promoters(txdb, upstream=2000, downstream=400)
PR

EX <- exons(txdb)
EX[1:4]
length(EX)
length(GR)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


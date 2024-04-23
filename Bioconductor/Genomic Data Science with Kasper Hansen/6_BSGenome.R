library(BiocManager)
library(BSgenome)
library(Biostrings)
library(IRanges)
library(GenomicRanges)
library(dplyr)
library(tidyverse)
library(plyranges)
library(ShortRead)
library(Rqc)
library(AnnotationHub)

# list of available BSgenome datasets---------------------------------------------------
available.genomes()

# List of Installed BSgenomes data packages--------------------------------------------
installed.genomes()

# Formal Definition of the Class--------------------------------------------------------------
showClass("BSgenome")

# Summary of Accessor-------------------------------------------------------------------------
.S4methods(class = "BSgenome")
.S4methods(class = "GenomeDescription")
# For Subclasses
showMethods(classes = "GenomeDescription", where = search())
showMethods(classes = "BSgenome", where = search())

# Get Help for BSgenome datasets-------------------------------------------------------------
help(BSgenome)

## Install-BSGenome-resoureces---------------------------------------------------------------
# Install the BSgenome package data (BSgenome.Scerevisiae.UCSC.sacCer2)  
# To install this package, start R (version "4.1")
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer2")
BiocManager::install("BSgenome.Celegans.UCSC.ce2")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Scerevisiae.UCSC.sacCer2)
library(BSgenome.Celegans.UCSC.ce2)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Celegans.UCSC.ce2
yeast <- BSgenome.Scerevisiae.UCSC.sacCer2
hg38 <- BSgenome.Hsapiens.UCSC.hg38
yeast
genome
hg38
Scerevisiae
metadata(hg38)
organism(hg38)
provider(hg38)

# What is yeast's main class?------------------------------------------------------
class(yeast)  # "BSgenome"
class(hg38)

# What is yeast's other classes?--------------------------------------------------
is(yeast)  # "BSgenome", "GenomeDescription"
is(hg38)

##--------------S4-checking----------------------------------
#------------Is yeast an S4 representation?--------------------------------
isS4(yeast)  # TRUE
str(yeast)

#---------------Seqname-Chromosome-Names-------------------------------------------
names(yeast)
seqnames(yeast)
names(hg38)
seqnames(hg38)

#------------object-summary---------------------------------------------
show(genome)
show(hg38)

#---------------------------seqinfo-Get-Sequence-Informations-------------------------
seqinfo(yeast)
yeast$chrI
yeast$chrXVI

# The BSgenome data packages only provide the forward strand sequence of every chromosome. 
# See ?reverseComplement for more details about the reverse complement of a DNAString object. It is important to remember that, in practice, the reverse strand sequence is almost never needed.The reason is that, in fact, a reverse strand analysis can (and should) always be transposed
# into a forward strand analysis
?reverseComplement

#---------------------------Chromosome-Number-in-genome--------------------------------
length(yeast)

#------------------------Sequence-Lengths-of-each-chromosome---------------------------
seqlengths(yeast)

#--------------Get Sequence Informations-------------------------------------
seqinfo(yeast)

#------------Get Chromosome M into an Object chrM and do some analysis--------------------
chrM <- yeast$chrM
length(chrM)
nchar(chrM)
sample(chrM)
reverse(chrM)
reverseComplement(chrM)
subseq(chrM, start = 85769, width = 10)

#------------Find Alphabet-&-Alphabet-Frequency for Chromososme M--------------
chrI <- genome$chrI
af <- alphabet(chrI)
af
af1 <- alphabetFrequency(chrI)
af1
sum(af1) == length(chrI)

#--------------LetterFrequency---------------------------------------------
letterFrequency(Scerevisiae$chrI, "CG")
letterFrequency(Scerevisiae$chrI, "CG", as.prob = TRUE)

##-----------------getseq----------------------------------------------------------------
#Genomes are often big, but interest usually lies in specific regions of them. Therefore, we need to subset a genome by extracting parts of it. To pick a sequence interval, use getSeq() and specify the name of the chromosome and the start and end of the sequence interval.
#S4 method getSeq() requires a BSgenome object
getSeq(yeast)
# OR,
# Display Chromosome with the following way- 
yeast$chrXV
# Select Chromosome sequences by name, one or many
getSeq(yeast, "chrI")
# Select start, end and or width 
# end = 10, select first 10 base pairs of each chromosome
getSeq(yeast, end = 10)
# Get the first 30 bases of chrM
getSeq(yeast, names = "chrM", end = 30)
getSeq(genome, "chrI", end = 10)
getSeq(yeast, names = "chrI", start = 100, end = 150)
getSeq(genome, "chrM", start = 100, end = 150)


##---------------Pattern-Coutn-Matching-----------------------------------------------
#------------------Count all exact matched of Pattern "ACCCAGGGAC"-------------
chrM <- genome$chrM
p1 <- "TTTTTT"
countPattern(p1, chrM)

#Like most pattern matching functions in Biostrings, the countPattern and matchPattern functions support inexact matching. One form of inexact matching is to allow a few mismatching letters per match. Here we allow at most one:
countPattern(p1, chrM, max.mismatch=1)
m1 <- matchPattern(p1, chrM, max.mismatch = 1)
m1
class(m1)
mismatch(p1, m1[4:5])

#-----------------pattern-matching---------------------------------------------
chrI <- genome$chrI
p2 <- DNAString("AAGCCTAAGCCTAAGCCTAA")
m2 <- matchPattern(p2, chrI)
m2
m2 <- matchPattern(p2, chrI, max.mismatch=2)
head(m2)
m2[1:4]
p2 == m2[1:4]
mismatch(p2, m2[1:4])
#list of exact matches and the list of inexact matches can both be obtained with:
m2[p2==m2]
m2[p2!=m2]
#Note that the length of m2[p2 == m2] should be equal to countPattern(p2, chrI, max.mismatch=0)
countPattern(p2, chrI)
countPattern(p2, chrI, max.mismatch = 2)
#Match a single string agains a set of sequence
vmatchPattern(p2, genome)
vmatchPattern(p2, genome)

## ----reverseComplement--------------------------------------------------------
p2 == reverseComplement(p2)

#Some views functionality of genome
ranges(m2)
genome$chrI[5:24]
genome$chrI[14445787: 14445806]
alphabetFrequency(m2)
genome$chrI[ start(m2):end(m2) ]
#Shift the m2 string to 10 bases ahed
shift(m2, 10)

## ----viewsVMatchPattern--------------------------------------------------
gr <- vmatchPattern(p2, genome)
gr
vi2 <- Views(genome, gr)
vi2

## ----annotationHub-------------------------------------------------------
ahub <- AnnotationHub()
qh <- query(ahub, c("sacCer2", "genes"))
qh
genes <- qh[[which(qh$title == "SGD Genes")]]
genes

## ----promoterGCcontent---------------------------------------------------
prom <- promoters(genes)
prom
head(prom, n = 3)

## ----promoterGCcontent2--------------------------------------------------
prom <- trim(prom)
promViews <- Views(Scerevisiae, prom)
promViews
gcProm <- letterFrequency(promViews, "GC", as.prob = TRUE)
head(gcProm)

## ----genomeGC-content-----------------------------------------------------------
params <- new("BSParams", X = Scerevisiae, FUN = letterFrequency, simplify = TRUE)
gccontent <- bsapply(params, letters = "GC")
gcPercentage <- sum(gccontent) / sum(seqlengths(Scerevisiae))
gcPercentage

## ----plotGC, fig=TRUE----------------------------------------------------
plot(density(gcProm))
abline(v = gcPercentage, col = "red")

#--------find-palindromic-sequence--------------------------------------
chrI <- genome$chrI
findPalindromes(chrI)

##--------------Compute-GC-Content-of-Entire-Genome--------------------------------
#--------BSapply-------------------------------
param <- new("BSParams", X = Scerevisiae, FUN = letterFrequency)
bsapply(param, letters = "GC")
head(bsapply(param, letters = "GC"))
unlist(bsapply(param, letters = "GC"))
#Calculate GC for entire genome
sum(unlist(bsapply(param, letters = "GC"))) / sum(seqlengths(Scerevisiae))
#GC content for individual chromosome
unlist(bsapply(param, letters = "GC", as.prob = TRUE))










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
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer2")
library(BSgenome.Scerevisiae.UCSC.sacCer2)
yeast <- BSgenome.Scerevisiae.UCSC.sacCer2

#Package Name: Biostring
#Description: Efficient manipulation of biological strings
#Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/Biostrings.html
#Documentation: PDF

#=============================Biostring==============================================
# Introduction to Biostrings
BiocManager::install("Biostrings")
??Biostrings

# Biological String Container
showClass("XString")
showClass("XStringSet")
showClass("BString")
showClass("BStringSet")

## ----IUPAC---------------------------------------------------------------
#The International Union of Pure and Applied Chemistry (IUPAC) codes: http://genome.ucsc.edu/goldenPath/help/iupac.html
IUPAC_CODE_MAP
# Biostrings Alphabet
# DNA 4 Bases
DNA_BASES
# RNA 4 Bases
RNA_BASES
# 20 Amino Acids
AA_STANDARD
# IUPAC_Code_Map
DNA_ALPHABET
RNA_ALPHABET
AA_ALPHABET

## ----DNAString, error=TRUE-----------------------------------------------
dna1 <- DNAString("ACGT-G")
dna1
DNAStringSet("ADE")
dna2 <- DNAStringSet(c("ACGT", "GTCA", "GCTA"))
dna2
#DNAString Subset---------------------------------
dna1[1:2]
dna2[1:2]
dna2[[2]]

## ----DNAStringSetNames---------------------------------------------------
names(dna2) <- paste0("seq", 1:3)
dna2

## ----Basic-Function-----------------------------------------------------------
width(dna2)
sort(dna2)
rev(dna2)
rev(dna1)
reverse(dna2)
complement(dna1)
reverseComplement(dna1)
complement(dna2)
dna2 <- unlist(dna2)

## ----counting------------------------------------------------------------
alphabetFrequency(dna1)
alphabetFrequency(dna2)
letterFrequency(dna2, letters = "GC")
consensusMatrix(dna2, as.prob = TRUE)
dinucleotideFrequency(dna2)
trinucleotideFrequency(dna2)

#Transcription Process (DNA to RNA)
dna_seq <- DNAString("ATGATCTCGTAA")
dna_seq
#DNAString Subset
dna_seq[2:4]
complement(dna_seq)
rna_seq <- RNAString(dna_seq)
rna_seq

# Translation RNA to Amino Acids
RNA_GENETIC_CODE
rna_seq
aa_seq <- translate(rna_seq)
aa_seq

# Shortcut Translation DNA to Amino Acids
translate(dna_seq)

#------------Find Alphabet-&-Alphabet-Frequency for Chromososme M--------------
chrI <- yeast$chrI
af <- alphabet(chrI)
af
af1 <- alphabetFrequency(chrI)
af1
sum(af1) == length(chrI)

#power of Biostrings comes to light when manipulating much larger sequences.
# Biostrings: Subsetting or Manipulating the sequences with subseq and unlist command
# Unlist the set, select the first 100 seq, and assign to dna
dnaa <- subseq(chrM, start = 1, width = 100)
dnaa
dna <- unlist(subseq(chrM, start = 1, width = 100))
dna
dna_seqs <- subseq(unlist(chrM), end = 100)
dna_seqs
rna <- RNAString(dna)
rna






#=====================================================================================
## DataCamp Sequence Handling. From a single sequence to a set
# Putting chrI chromosome into a vector
chrI <- yeast$chrI

# Create a new set from a single sequence
yeastset <- DNAStringSet(chrI, start = c(1, 101, 201, 301), end = c(100, 200, 300, 400))
yeastset
length(yeastset)
width(yeastset) # Width only used on set and gives number of character per sequenc

# Collate all the sequences in a single sequence using unlist() function, and for single seq width will not work here.
yeastt <- unlist(yeastset)
length(yeastt)
nchar(yeastt)

# Sort a sequence set, (if decreasing = TRUE - sorting Z to A order), (if decreasing = FALSE - sorting A to Z order )
sort(yeastset, decreasing = FALSE)

##-------------------Sequence-Manipulations-----------------------------------
# How to create complementary sequence
a_seq <- DNAString("ATCTGATGCCGG")
complement(a_seq)
# How to reverse the sequence set or change the order of sequence
yeastset
rev(yeastset)
# Reverse the each sequence or change the order of sequence in a set from right to left.
reverse(yeastset)
# Reverse complement
reverseComplement(yeastt) 
#OR, 
reverse(complement(yeastt))

#===============================================================================================
# Pattern matching in Sequences
# Find pattern is a set of sequences using vmatchPattern or vcountPattern
vmatchPattern(pattern = "CCCA", subject = yeastset, 
              max.mismatch = 0)

# Find pattern in a single sequence
matchPattern(pattern = "CCCA", 
             subject = yeastt, max.mismatch = 1)

# Find the palindromes in yeastt
findPalindromes(yeastt)


## PDF - Count all exact matched of Pattern "ACCCAGGGAC"
getSeq(yeast, "chrM")
p1 <- "TTCATAATT"
vcountPattern(p1, chrM)

# Like most pattern matching functions in Biostrings, the countPattern and matchPattern functions support inexact matching. One form of inexact matching is to allow a few mismatching letters per match. Here we allow at most one:
vcountPattern(p1, chrM, max.mismatch=1)
m1 <- vmatchPattern(p1, chrM, max.mismatch = 1)
m1
class(m1)
mismatch(p1, m1[4:5])




#--------Manipulating collections of GRanges-----------------------------------------
#---------slidingwindows-----------------
# SlidingWindows are useful to split GRanges object into sub-elements
# "width" is the total number of letters for each new ranges
# "step" is the distance between the ranges
# Genes has been split into new ranges of widht 20000 bases and distance between the ranges is 10000 because of step
# Each range has 10000 bases overlap because of "width-step". The last range could be shorter. 
slidingWindows(hg_chrXG, width = 20000, step = 10000)

# Prefilter chromosome X by seqlevels()
seqlevels(hg) <- c("chrX")
# Transcripts Extraction
transcripts(hg, columns = c("tx_id", "tx_name"), filter = NULL)
# Exons Extraction 
exons(hg, columns = c("tx_id", "tx_name"), filter = list(tx_id = "228616"))
# Exons by transcript
exonsBytx <- exonsBy(hg, by = "tx")
exonsBytx
# Trasncript ID
abcd1_179161 <- exonsBytx[["235985"]]
abcd1_179161
width(abcd_179161)














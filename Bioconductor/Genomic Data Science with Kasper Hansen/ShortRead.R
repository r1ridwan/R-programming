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

#Package Name: ShortRead
#Description: FASTQ input and manipulation
#Bioconductor URL: https://bioconductor.org/packages/release/bioc/html/ShortRead.html
#Documentation URL: PDF

BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer2")
library(BSgenome.Scerevisiae.UCSC.sacCer2)
yeast <- BSgenome.Scerevisiae.UCSC.sacCer2

#=============ShortRead=============================================================
## Introduction to ShortReads (FASTA, FASTQ, readFasta, writeFasta, readFastq, writeFastq,)
## the main functionalities of this package are:Data input, Quality assessment, Data transformation, Access to downstream analysis
BiocManager::install("ShortRead", force = TRUE)
help("ShortRead")

# Method Accessor for FASTQ and FASTA format
methods(class = "ShortRead")
methods(class = "ShortReadQ")

#------------Read-Write-FASTA-FASTQ-File-------------------------------------------------
# Reading FASTA file S4 method-1 in fasta folder
fasample <- readFasta(dirPath = "fasta/", pattern = "fasta")
fasample
class(fasample)
fastqq <- readFastq(dirPath = "fasta/", pattern = "SRR1971253_pass.fastq")
fastqq
class(fastqq)
sread(fastqq)[1:2]
id(fastqq)[1:2]
quality(fastqq)[1:2]

# Reading FASTA Files in fastafolder method-2
dna <- readDNAStringSet("fasta/sequence.fasta")

#Write FASTA file
writeFasta(fasample, file = "fasta/fastaaa.fasta")
#Write FASTQ file
writeFastq(fastqq, file = "fasta/fastqqqfile.fastq")

# Read FASTQ File in S4 method and check other functions (SRR1971253 SRA file was downloaded for analysis)
# Arabidopsis Thaliana plant SRA file used here
fasttq <- readFastq(dirPath = "fasta/", pattern = "fastq")
fasttq
class(fasttq)
class(sread(fasttq))
id (fasttq)
# Check widths of each read
rds <- sread(fasttq)
wid <- rds@ranges@width
wid

# Read FASTQ file 
f = readFastq("SRR1971253_pass.fastq")
head(f)
summary(f)
reads <- sread(f)
head(reads)
# Check width of each read
widths <- reads@ranges@width
widths 


##---------Set.Seed------------------------------------------------------
# Set the seed to draw the same read sequence every time
set.seed(1234)

# Use FastqSampler with f(file path) and select 100 reads from 60836 reads and each read with 50bp
# Extract a subsample from fastq file
fs <- FastqSampler("fasta/SRR1971253_pass.fastq", 100)
fs
fs2 <- FastqSampler(con = "fasta/SRR1971253_pass.fastq", n = 200) #con = "file path", n= number of reads
fs2
# Generate new sample yield
#Save the yield of 100 read sequenc
my_sample <- yield(fs)
my_sample
head(my_sample)
class(my_sample)
length(my_sample)


# Count FASTQ file in S4 method
ct = countFastq(dirPath = "SRR1971253_pass.fastq", pattern = "fastq")
ct

#============Sequence Quality====================================================
## Sequence Quality
# Endcoding sequence Quality or FASTQ sequence quality score
encoding(quality(fasttq))
quality(fasttq)

sread(fasttq)[1]
# Quality is represented as ASCII Characters
quality(fasttq)[1]
# Phred Quality Instance
qp <- PhredQuality(quality(fasttq))
# Transform encoding phredquality into scores (Score 30 is considered as good quality score)
qs <- as(qp, "IntegerList")
qs
#Convert Quality
as(quality(fastqq), "matrix")[1:2, 1:10]

##------------Sequence-Quality-Assessment----------------------------------------
# Quality Assessment 
qaSummary <- qa(fasttq, type = "fastq", lane = 1)
qaSummary
names(qaSummary)

#QA elements can be accessed with qaSummary[["nameElement"]], where nameElement is the name of the element you wish to inspect.
qaSummary[["readCounts"]]
qaSummary[["alignQuality"]]
qaSummary[["sequenceDistribution"]]

#Quality Score Encoding by Illumina
#URL: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm

# Get HTML Report
browseURL(report(qaSummary))

#Check Nucleotide frequency per cycle
# Check ShortRead Sequence Alphabet
alphabet(sread(fasttq))
abc <- alphabetByCycle(sread(fasttq))
abc

# Each observation is a letter and each vaiable is a cycle. First, Select the first 4 rows nucleotides A, C, G, T
# Then Transpose
nucbycycle <- t(abc[1:4, ])
nucbycycle
# Convert to tibble and add cycle numbers
nucbycycle <- nucbycycle %>% as.tibble() %>%
  mutate(cycle = 1:50)
nucbycycle

#------------Frequency-Plotting--------------------------------------------------
#Nucleotide Frequency Plot
# Create a line plot of cycle vs. count
glimpse(nucbycycle)
nucbycycle %>% 
  # Gather the nucleotide letters in alphabet and get a new count column
  gather(key = alphabet, value = count , -cycle) %>% 
  ggplot(aes(x = cycle, y =  count, color = alphabet)) +
  geom_line(size = 0.5 ) +
  labs(y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

#=========Match-and-Filter=======================================================
## Matching and Filtering 
# Counting Duplicates- TRUE is the number of duplicates
table(srduplicated(fasttq))

# Cleaning reads from duplicates
cleanRead <- fasttq [srduplicated(fasttq) == FALSE]
# Counting Duplicates
table(srduplicated(cleanRead))

# Use a custom filter to remove reads from sequence
# This filter to remove reads shorter than a minimum number of bases
minWidth <- 51
readWidthCutOff <- srFilter(function(x) {width(x) >= minWidth},name = "MinWidth")
fasttq[readWidthCutOff(fasttq)]

# only in those reads that start with the pattern "ATGCA". 
#A tiny filtering function can do the job, making use of the srFilter() function:
myStartFilter <- srFilter(function(x) substr(sread(x), 1, 5) == "ATGCA")
fasttq[myStartFilter(fasttq)]

# Check class of fqsample
class(fasttq)
# Filter reads into selectedReads using myStartFilter
selectedReads <- fasttq[myStartFilter(fasttq)]
# Check class of selectedReads
class(selectedReads)

# Save your filter
myfilter <- nFilter(threshold = 10, .name = "cleaNFilter")
# Use the filter at reading point
filtered <- readFastq(dirPath = "fasta/",
                      pattern = ".fastq",
                      filter(myfilter))
# You will retrieve only those reads that have a maximum of 10 Nucleotides
filtered

#ID filter example
myidfilter <- idFilter(regex = ":3:1") #It will return only those ids contain only regular expression
#Optional parameters are .name, fixed and exclude
#Use the filter at reading point
filtered <- readFastq(dirPath = "fasta/", pattern = ".fastq",
                      filter = myfilter)

filtered

#Filter to remove poly-A regions
# Will return the sequence that have a maximum number of 10 consecutive A's
myfilterPolyA <- polynFilter(threshold = 10, nuc = c("A")) 
filtered[myfilterPolyA(filtered)]

















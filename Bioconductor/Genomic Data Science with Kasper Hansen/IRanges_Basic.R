##===============Bioconductor for genomic data science================================
#URL: http://kasperdanielhansen.github.io/genbioconductor/html/IRanges_Basic.html
## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(IRanges)
library(plyranges)

## ---------------------IRanges------------------------------------------------------------
# IRanges construction with different start, end and width combinations
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir1
ir2 <- IRanges(start = c(1,3,5), width = 3)
all.equal(ir1, ir2)
ir1 == ir2
ir3 <- IRanges(start = 1:10, width = 10:1)
ir3
ir4 <- IRanges(start = 1:10, end = 11)
ir4
ir5 <- IRanges(start = 1:10, width = 10)
identical(ir3, ir4) && identical(ir3, ir5)
ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40),
              width = c(12, 6, 6, 15, 6, 2, 7))
ir
## IRanges and Genomic Structures
# IRanges with numeric arguments
myIranges <- IRanges(start = 20, end = 30)
myIranges
myIRanges_width <- IRanges(start = c(1, 20), width = c(30, 11))
myIRanges_width
(myIranges_end <- IRanges(start = c(1, 20), end = 30))

#---------------IRanges-with-Names-------------------------------------
seq_2 <- IRanges(start = c(5, 35, 50),
                 end = c(12, 39, 61),
                 names = LETTERS[1:3])
seq_2

## ----ir_width------------------------------------------------------------
start(ir1)
end(ir1)
width(ir)
width(ir1) <- 5
ir1

#----------------------------Rle---------------------------------------
#Rle stands for Run length encoding
#Computes and stores the lengths and values of a vector or factor
#Rle is a general S4 container use to save long repetitive vectors efficiently
some_numbers <- c(3, 2, 2, 2, 3, 3, 4, 2)
Rle(some_numbers)
#IRanges with logical Rle
gi <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE)
myRle <- Rle(gi)
myRle
IRanges(start = myRle)#True element in myRle
# IRnum1: start - vector 1 through 5, end - 100  
IRnum1 <- IRanges(start = c(1, 2, 3, 4, 5), end = 100)
IRnum2 <- IRanges(start = c(11, 90), width = c(89, 10))
IRlog1 <- IRanges(start = Rle(c(F, T, T, T, F, T, T, T)))
# Print objects in a list
print(list(IRnum1 = IRnum1, IRnum2 = IRnum2, IRlog1 = IRlog1))

##--------------Subsetting an IRanges object----------------------- 
ir[1:4]
ir[start(ir) <= 15]

## ----ir_names------------------------------------------------------------
names(ir) <- paste("A", 1:7, sep = "")
ir

## ----ir_dim--------------------------------------------------------------
dim(ir1)
length(ir1)

## ----ir_subset-----------------------------------------------------------
ir[1]
ir["A1"]

## ----------concatenate and giving names---------------------------------------------------------
iname = c(ir1, ir2)
names(iname) <- paste("A", 1:6, sep = "")
iname

## ----irNormal1, echo=FALSE-----------------------------------------------
ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
ir
## ----irNormal2, echo=FALSE, fig.height=2, small.mar=TRUE-----------------
plotRanges(ir)

## ----irNormal3, echo=FALSE, fig.height=1.75, small.mar=TRUE--------------
plotRanges(reduce(ir))

## ----irNormal4-----------------------------------------------------------
ir
reduce(ir)

## ----irDisjoin1, eval=FALSE----------------------------------------------
## disjoin(ir1)

## ----irDisjoin2, echo=FALSE, fig.height=2, small.mar=TRUE----------------
plotRanges(ir)

## ----irDisjoin3, echo=FALSE, fig.height=1.75, small.mar=TRUE-------------
plotRanges(disjoin(ir))

## -----------ir_resize----------------------------------------------------
ir
resize(ir, width = 1, fix = "start")
resize(ir, width = 1, fix = "center")
resize(ir, width = 5, fix = "end")

## ----ir_sets_union_and_intersect------------------------------------------
ir1 <- IRanges(start = c(1, 3, 5), width = 1)
ir1
ir2 <- IRanges(start = c(4, 5, 6), width = 1)
ir2
union(ir1, ir2) #or
reduce(c(ir1, ir2))
intersect(ir1, ir2)
ir3 <- IRanges(start = 1:10, width = 11)
ir3
ir4 <- IRanges(start = 11:20, width = 11)
ir4
union(ir3, ir4)
intersect(ir3, ir4)

## ----union2--------------------------------------------------------------
reduce(c(ir1, ir2))

## ----findOverlaps--------------------------------------------------------
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width = 3)
ov <- findOverlaps(ir1, ir2)
ov
ov2 <- findOverlaps(ir3, ir4)
ov2
ov3 <- findOverlaps(ir3, reduce(ir3))
ov3
as.matrix(ov2)
as.matrix(ov3)
?findOverlaps

## ----findOverlaps_ill----------------------------------------------------
intersect(ir1[subjectHits(ov)[1]],
          ir2[queryHits(ov)[2]])

## ----subjectHits---------------------------------------------------------
queryHits(ov2)
unique(queryHits(ov2))

## ----argsFindOverlaps, tidy=TRUE-----------------------------------------
args(findOverlaps)

## ----countOverlaps-------------------------------------------------------
countOverlaps(ir1, ir2)
countOverlaps(ir3, ir4)
cov <- coverage(ir1)
cov
cov2 <- as.vector(cov)
cov2

## ----nearest-------------------------------------------------------------
ir1
ir2
nearest(ir1, ir2)
nearest(ir3,ir4)

#-------------Adjusting strts, end and width--------------------------------
ir_new <- IRanges(start = c(1:3),end=c(5,2,8))
ir_new
# Shift start and end value with the following command
shift(ir_new, 10)
narrow(ir3, start = 1:5, width = 5)

##----------------------Plyranges------------------------------------------------
# Modify (or Mutate) a genomic intervals by altering the width
# Use the mutate verb with anchor_*() adverbs
# To work the mutate function 'plyranges' package is requrired
ir_new <- IRanges(start = c(1:3), end = c(5, 2, 8))
ir_new
mutate(ir_new, width = 10) # start value thik rekhe end value change hobe and width 10 hobe
mutate(anchor_start(ir_new), width = 10) # ei function a start value thik rekhe end value change hoise jate width 10 hoy
mutate(anchor_end(ir_new), width = 10) # ei function a end value thik rekhe start value change hoise jate kore width 10 hoy
## Stretch a genomic interval by an intger amount 
ir_new
stretch(anchor_center(ir_new), extend = 10) # ekhane total width 10 barbe, start 5 and end 5 barbe
# Shifting Genomic regions or bases in left direction or right direction
shift_left(ir_new, shift = 10)
shift_right(ir_new, shift = 10)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


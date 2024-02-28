##========Coursera====================================================
##============Bioconductor for Genomic Data Science==================
## ----dependencies, warning=FALSE, message=FALSE--------------------------
## URL: http://kasperdanielhansen.github.io/genbioconductor/html/R_Base_Types.html

## ----numeric-------------------------------------------------------------
x <- 1:10
names(x) <- LETTERS[1:10]
class(x)
x[1:10]
x[c("A", "B")]
class(x)


## ----uNames--------------------------------------------------------------
x <- 1:10
names(x) <- c("A", "A", "B", "C", "E", "D", "C", "A", "B", "D")
x
x["A"]

## ----uNames2-------------------------------------------------------------
anyDuplicated(names(x))
names(x) <- c("A", "B", "C")
anyDuplicated(names(x))

## ----intNum--------------------------------------------------------------
x <- 1
class(x)
x <- 1:3
class(x)

## ----intNum2-------------------------------------------------------------
x <- 1L
class(x)
as.numeric(x)

## ----machine-------------------------------------------------------------
.Machine$integer.max
2^31 -1 == .Machine$integer.max
round(.Machine$integer.max / 10^6, 1)
# This number is smaller than number of bases in the human genome. To fix this use as.numeric(). 

## ----matrices------------------------------------------------------------
x <- matrix(1:9, ncol = 3, nrow = 3)
rownames(x) <- c("A","B","C")
rownames(x) <- LETTERS[1:3]
colnames(x) <- letters[1:3]
x
dim(x) #dimension of matrix
nrow(x)
ncol(x)

##-----subsetting matrices-----------------------------------------------------------
#tow dimensional subsetting of matirces, [rows, columns]
x[1:2,]
x["B", ]
x[x >= 5]

## ----matrixSubset2-------------------------------------------------------
x[1,]
x["A", ] # it gives the output as a vector not as a matrix
x[1,,drop = FALSE] # it gives the output as a matrix


## ----createMatrix--------------------------------------------------------
matrix(1:9, ncol = 3, nrow = 3)
matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)

## ----list----------------------------------------------------------------
x <- list(1:3, letters[1:3], is.numeric)
x
names(x) <- c("numbers", "letters", "function")
x[1:2]
x[1]
x[[1]]

## ----list2---------------------------------------------------------------
x$letters
x["letters"]
x$let
x["let"]

## ----as.list-------------------------------------------------------------
as.list(1:3)
list(1:3)

## ----lapply--------------------------------------------------------------
x <- list(a = rnorm(3), b = rnorm(3))
lapply(x, mean)

## ----sapply--------------------------------------------------------------
sapply(x, mean)

## ----df------------------------------------------------------------------
x <- data.frame(sex = c("M", "M", "F"), age = c(32,34,29))
x

## ----dfColumns-----------------------------------------------------------
x$sex
x[["sex"]]

## ----df2-----------------------------------------------------------------
x <- data.frame(sex = c("M", "M", "F"), age = c(32,34,29), stringsAsFactors = FALSE)
x$sex

## ----dfApply-------------------------------------------------------------
sapply(x, class)

## ----as.X, error=TRUE----------------------------------------------------
x
as.matrix(x)
as.list(x)

## ----as------------------------------------------------------------------
library(methods)
as(x, "matrix")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()


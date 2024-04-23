## Course: Bioconductor for Genomic Data Science by Kasper Hanssen
## Hansen Lab: https://www.hansenlab.org/
## Full Course: https://kasperdanielhansen.github.io/genbioconductor/
## HTML URL: http://kasperdanielhansen.github.io/genbioconductor/html/R_Base_Types.html
## YouTube: https://www.youtube.com/watch?v=x75azcLzQzs


## R - Base Types
# 1. numeric-----------------
x <- 1:10
names(x) <- LETTERS[1:10]
class(x)
x[1:10]
x[c("A", "B")]
class(x)

# uNames----------------
# All vectors can have missing values.
# Note: names of vectors does not need to be unique. This can lead to subsetting problems:
x <- 1:10
names(x) <- c("A", "A", "B", "C", "E", "D", "C", "A", "B", "D")
x
x["A"]

# uNames2----------------
# check for unique names
anyDuplicated(names(x))
names(x) <- c("A", "B", "C")
anyDuplicated(names(x))
duplicated(names(x))
unique(names(x))


# 2. integerNum-----------------
x <- 1
class(x)
x <- 1:3
class(x)

# intNum2----------------
x <- 1L
class(x)
as.numeric(x)

# machine-------------
# The maximum integer
.Machine$integer.max
2^31 -1 == .Machine$integer.max
round(.Machine$integer.max / 10^6, 1)
# This number is smaller than number of bases in the human genome. To fix this use as.numeric(). 


# 3. matrices----------------
# matrix is a two-dimensional object
x <- matrix(1:9, ncol = 3, nrow = 3)
rownames(x) <- c("A","B","C")
rownames(x) <- LETTERS[1:3]
colnames(x) <- letters[1:3]
x
dim(x) #dimension of matrix
nrow(x)
ncol(x)

# create Matrix----------------
matrix(1:9, ncol = 3, nrow = 3)
matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)

# sub-setting matrices--------------
# two dimensional sub-setting of matrices
# the first dimension is rows and the second is columns
x[1:2,]
x["B", ]
x[x >= 5]

# matrix Subset2---------------
x[1,]
x["A", ] # it gives the output as a vector not as a matrix
x[1,,drop = FALSE] # it gives the output as a matrix


# 4. list-------------
x <- list(1:3, letters[1:3], is.numeric)
x
names(x) <- c("numbers", "letters", "function")
x[1:2]
x[1]
x[[1]]

# list2---------------
x$letters
x["letters"]
x$let
x["let"]

# as.list------------
as.list(1:3)
list(1:3)

# lapply--------------
x <- list(a = rnorm(3), b = rnorm(3))
lapply(x, mean)

# sapply----------------
sapply(x, mean)


# 5. Data frames---------------
x <- data.frame(sex = c("M", "M", "F"), age = c(32,34,29))
x

# dfColumns-----------
x$sex
x[["sex"]]

# df2--------------
x <- data.frame(sex = c("M", "M", "F"), age = c(32,34,29), stringsAsFactors = FALSE)
x$sex

# dfApply---------------
sapply(x, class)


# 6. Conversion-------------
# as.X, error=TRUE
as.matrix(x)
as.list(x)

# as
library(methods)
as(x, "matrix")


# sessionInfo, echo=FALSE-----------
sessionInfo()



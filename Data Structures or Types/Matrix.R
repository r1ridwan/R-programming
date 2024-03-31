## Matrix in R
# a matrix is a two-dimensional array-like structure
matrix(1:9)
matrix(1:9, ncol = 3)
matrix(1:9, ncol = 3, nrow = 3)
matrix(1:9, ncol = 3, nrow = 3, byrow = T)

mat <- matrix(1:9, ncol = 3, nrow = 3, byrow = T)
mat

# subset matrix
mat[1]

# subset entire row
mat[1,]

# subset column
mat[,1]

# Matrix with dimension names
mat2 <- matrix(1:12, ncol = 3, nrow = 4, byrow = T, 
              dimnames = list(c("row1", "row2", "row3", "row4"),
                              c("col1", "col2", "col3")))

## Matrix Properties
# Get column names
colnames(mat2)
# get row names
rownames(mat2)
# dimensions
dim(mat2)


## Matrix addition
mat3 <- matrix(9:1, nrow = 3, ncol = 3, byrow = T)
result <- mat + mat3

# Display the result of addition
print(result)


# Transpose of the matrix
t(mat)








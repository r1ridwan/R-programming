# 1. Information about variables
x <- 6
n <- 1:4
let <- LETTERS[1:4]
df = data.frame(n, let)

    
# 2. Information About Existence
# List currently defined variables
ls()

# Check if a variable named "x" exists
exists("x")

# Check if "y" exists
exists("y")

  
# 3. Information about size/structure
# Get information about structure
str(n)
str(df)

# Get the length of a vector
length(n)

# Length probably doesn't give us what we want here:
length(df)

# Number of rows
nrow(df)

# Number of columns
ncol(df)

# Get rows and columns / dimensions
dim(df)


# 4. How to check and extract a variable from a data set
# install xlsx package
install.packages("xlsx")
library(xlsx)

# Loading LungCapData file
LungCapData = read.xlsx(file="LungCapData.xlsx", sheetIndex = 1)

# Check structure of dataset
str(LungCapData)

# Getting names of this dataset 
names(LungCapData)

# Getting mean for Age variable ~ it shows errors before we attach it
mean(Age)

# We have to use dollar sign to detect the variable on this dataset
mean(LungCapData$Age)

# Lets attach the data in R's memory
attach(LungCapData)

# Now we don't need to use dollar sign to detect variables
mean(Age)
mean(Height)

# We can detach the data also
detach(LungCapData)

# Now it will show error again
mean(Age)


# 5. Getting the type of variable with class command
# Attach the data again
attach(LungCapData)
names(LungCapData)
class(LungCap)
class(Age)
class(Height)
class(Smoke)
class(Gender)
class(Caesarean)

# To find out what levels or category have on this character variable

# get the summary 
summary(LungCapData)

  
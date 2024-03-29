library('tidyverse')

# Functions
# for identifying duplicated elements and
duplicated(): 
# for extracting unique elements
unique():
# [dplyr package] ~ to remove duplicate rows in a data frame.
distinct() 


# Remove Duplicates
# Given the following vector:
x <- c(1, 1, 4, 5, 4, 6)
# To find the position of duplicate elements in x, use this:
duplicated(x)
# Extract duplicate elements:
x[duplicated(x)]
# remove duplicated elements
x[!duplicated(x)]


# Remove duplicates rows from a data frame based on a column values
# Load data
data <- as.tibble(iris)
duplicated(iris)
# Find the duplicated row number
which(duplicated(iris)) 
# Remove duplicates based on Sepal.Width columns
# ekhane 'speal.width' jei data gula duplicate hoise oi full row gula delete hoye jabe. 
data = data[!duplicated(data$Sepal.Width), ] 


# Extract unique elements
x = c(1,1,4,5,4,6)
unique(x)

unique(data)
data2 = unique(data)


# Remove duplicate rows in a data frame
# Remove duplicate rows based on all columns
data %>% distinct()

# Remove duplicate rows based on certain columns (variables):
# Remove duplicated rows based on Sepal.Length
data %>% distinct(Sepal.Length, .keep_all = TRUE)

# Remove duplicated rows based on Sepal.Length and Petal.Width
data3 = data %>% distinct(Sepal.Length, Petal.Length, .keep_all = TRUE)










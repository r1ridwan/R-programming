## Vector
# Vectors are the most basic R data objects and there are six types of atomic vectors. 
# logical, integer, double, complex, character and raw.

library(tidyverse)

typeof(letters)
typeof(1:12)

x <- list("a", "b", 1:10)
x
length(x)

# Vectors are generally created using the c() function
v <- c(1,4,4,3,2,2,3)
class(v)
v[c(2,3,4)]
v[1:3]

# To create vectors of consecutive numbers, the : operator is very helpful
v2 <- 1:20
class(v2)

# LinK: https://www.tutorialspoint.com/r/r_vectors.htm


# 1. Combining Vectors:
# http://www.r-tutor.com/r-introduction/vector/combining-vectors
# Vectors can be combined via the function c. For examples, the following two vectors n and s are combined into a new vector containing elements from both vectors.
n= c(2, 3, 5) 
s= c("aa", "bb", "cc", "dd", "ee") 
combined=c(s,  n)
print(combined)



# 2. SINGLE ELEMENT VECTOR Types
# Atomic vector of type CHARACTER
typeof("abc") 
v6 <- c("a", "b")
class(v6)

# Numeric vector
v4 <- c(0.5, 0.6)
class(v4)

# Atomic vector of type DOUBLE
typeof(12.5)

# Atomic vector of type INTEGER
typeof(63L)
v7 <- 1:10
class(v7)

# Atomic vector of type LOGICAL
typeof(TRUE)
v5 <- c(TRUE, FALSE)
class(v5)

# Atomic vector of type COMPLEX
typeof(2+4i)
v8 <- c(1+0i, 2+0i)
class(v8)

# Atomic vector of type RAW
typeof(charToRaw('hello'))

# Mixed vector
# Character
x1 <- c(1.7, "a")
class(x1)

# Numeric
y1 <- c(TRUE, 2)
class(y1)

# character
z1 <- c("a", F)
class(z1)



# 3. MULTIPLE ELEMENTS VECTOR
v1 <- 5:13
v1
v2 <- 3.8:11.4
v2



# 4. Using sequence (Seq.) operator
# Create vector with elements from 5 to 9 increment by 0.4
print(seq(5,9, by=0.4))
v3 <- seq(from = 1, to = 20, by = 3)
class(v3)
v3 <- seq(1, 20, 3)
v3

# Index sequence vector using seq_len()
index_seq <- seq_len(15)
index_seq



# 5. Repeated values vector using rep()
repeated_values <- rep(1, times = 5)
repeated_values



# 6. Vector created using vectorized operation
new_vector <- repeated_values * 2
new_vector
class(new_vector)
heights <- c(1.2, 1.3, 1.4, 1.5, 1.6) * 100
heights


## 7.  Sub-setting
# ACCESSING VECTOR ELEMENTS
ages <- c(22, 33, 45, 62, 34)
# subset/access by element position
ages[3]

# Accessing vector elements using multiple position
t <- c("Sun","Mon","Tue","Wed","Thurs","Fri","Sat")
u <- t[c(3,1,7)]
print(u)

# Accessing vector elements using logical indexing
v<- t[c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)] 
print(v)

# Accessing vector elements using 0/1 indexing.
y <- t[c(1,0,0,0,0,0,7)]
print(y)

# subset by sequence (: operator)
ages[1:3]

# reverse subset
ages[-1]

# subset by logical condition
vector <- c(10, 20, 30, 40, 50, 60)
vector[vector > 30]
vector[vector < 30]



# 8. VECTOR MANUPULATION
v1 <- c(3,8,4,5,0,11)
v2 <- c(4,11,0,8,1,2)

# Vector Addition 
add.result<- v1+v2 
add.result

# Vector subtraction 
sub.result<- v2-v1 
sub.result

# Vector Multiplication
multi.result<- v1*v2 
print(multi.result)

# Vector Division 
divi.result<-v1/v2 
print(divi.result)



# 9.  VECTOR RECYCLING
# If two vector are unequal in length, then shorter vector will recycle to match the longer vector, here below the shorter vector is v2, and it recycle two time to match the longer cycle. 
v1<- c(3, 4, 5, 6, 7, 8, 9)
v2<-c(4, 11)# v2 becomes (4,11,4,11,4,11,4)
add.result<- v1+v2
print(add.result)
sub.result<-v1-v2
print(sub.result)



# 10. VECTOR ELEMENTS SORTING
v<- c(3,8,4,5,0,11, -9, 304)
sort.result<- sort(v)
print(sort.result)
revsort.result<- sort(v, decreasing = TRUE)
print(revsort.result)
revsort.result<- sort(v, decreasing = FALSE)


v <- c("Blue", "Red", "Yellow", "violet", "Green")
sort.resutl <- sort(v)
sort.resutl
revsort.result<-sort(v, decreasing = FALSE)
print(revsort.result)
revsort.result<-sort(v, decreasing = TRUE)
print(revsort.result)



# 11. Vector Arithmetic 
# For example, suppose we have two vectors a and b.
a = c (2, 3, 4, 5, 7)
b = c (4, 7, 8, 9, 10)
# We can multiply a and b vector with 5
5 * a
5 * b
# We can sum this two vector together
a + b
# we can subtraction, multiplication and division of these two vector
a - b
a * b
a / b



# 12. Name Vector Members 
# We can assign names to vector members. For example, the following variable v is a character string vector with two members.
v <- c("Mary", "Nishi")
names(v)= c("First", "Last")
v["First"]
v["Last"]
v[c("Last", "First")]



# 13. Making a vector filled with values
# 1 will be repeated 50 times
rep(1, 50)
# 3 will be repeated 10 times
rep(3, 10)
# False will be output and for 50 times
rep(F, 50)
# 1 to 5 will be repeated for 4 times
rep(1:5, 4)

rep(1:5, each=4)
# Use it on a factor
rep(factor(LETTERS[1:3]), 5)



# 14. Vector Types Conversion or Data type conversion
# Treating strings as factors or characters.
# By default, strings  the data are converted to factors. If you load the data below with read.csv, then all the text columns will be treated as factors, even though it might make more sense to treat some of them as strings. If you don't want your data to be treated as factor instead of string then use the following command.
data <- read.csv("new_file.csv", stringsAsFactors = FALSE)

# You might have to convert some columns to factors
data$country <- factor(data$country)
class(data$country)



# Create a sample data frame
data <- read.table(header=T, text='
 subject sex size
       1   M    7
       2   F    6
       3   F    9
       4   M   11
 ')
data[1, 3]
data[3, 2]
data[1, "size"]

data[1:2, ]
data[c(1,2), ]

data[1:2, c(2,3)]
data[c(1,2), c("sex", "size")]
data[c(1,2), c(2,3)]


v <- c(1,4,4,3,2,2,3)
v > 2
v[v>2]
v[c(F,T,T,T,F,F,T)]

data <- read.table(header=T, text='
 subject sex size
       1   M    7
       2   F    6
       3   F    9
       4   M   11
 ')
data$subject <3
data[data$subject<3, ]
data[c(T, T, F, F), ]
# It is also possible to get the numeric indices of the TRUEs
which(data$subject < 3)

#Here is the vector again
v
#Drop the first element
v[-1]
# Drop the 4th element
v[-4]
#Drop the first three
v[-1:-3]
#Drop the just last element
v[-length(v)]
v = c (2, 3, 4, 5)
s = c ("two", "theree", "four", "five")
c(v,s)
c(s,v)


v1 <- c(3, 4, 5, 6, 7, 8, 9)
v2 <-c(4, 11)
add.result <- v2 + v1
print(add.result)
sub.result <- v2 - v1
print(sub.result)











# strings
str <- "Hello"
class(str)


str2 <- 'Hello'
class(str2)

# the following string will show error because of single quotation
str3 <- 'ridwan's

# its better to use double quotation in string
str4 <- "ridwan's"


# 1. string operations
# Find the length of a string
seq <- "AGGCTTGGGCCTAAA"
nchar(seq)

# join two string
s1 <- "ACCT"
s2 <- "GCCTA"
s <- paste(s1, s2)

# compare two strings
s1 == s2
s3 <- "ACCT"
s1 == s3


# 2. change the string case
name = "ridwan"

# uppercase
toupper(name)

# lowercase
tolower(name)












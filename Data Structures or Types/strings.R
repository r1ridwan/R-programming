# strings
str <- "Hello"
class(str)

str2 <- 'Hello'
class(str2)

# the following string will show error because of single quotation
str3 <- 'ridwan's

# its better to use double quotation in string
str4 <- "ridwan's"


# string operations
# 1. Find the length of a string
seq <- "AGGCTTGGGCCTAAA"
nchar(seq)


# 2. join two string
s1 <- "ACCT"
s2 <- "GCCTA"
s <- paste(s1, s2)


# 3. compare two strings
s1 == s2
s3 <- "ACCT"
s1 == s3


# 4. change the string case
name = "ridwan"

# uppercase
toupper(name)

# lowercase
tolower(name)












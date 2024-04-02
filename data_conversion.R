# 1. as.character(): Converts the input to character data type
# convert numeric to character
x <- c(12, 12, 134)
y <- as.character(x)
class(y)


# 2. as.numeric(): Converts the input to numeric data typ
# character to numeric
x <- c("10", "11", "12")
y <- as.numeric(x)
class(y)


# 3. as.integer(): Converts the input to integer data type
# numeric to integer
x <- 10.5
y <- as.integer(x)
class(y)


# 4. as.logical(): Converts the input to logical (boolean) data type
# numeric to logical
x <- 0
class(x)
y <- as.logical(x)
class(y)


# 5. as.factor(): Converts the input to a factor data type
x <- c("Male", "Female", "Male")
class(x)
y <- as.factor(x)
class(y)


# 6. as.list(): Converts the input to a list
x <- c(1, 2, 3)
class(x)
y <- as.list(x)
class(y)




















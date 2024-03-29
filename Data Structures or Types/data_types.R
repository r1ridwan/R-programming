# 1. Numeric
# Represents numeric values (integers or floating-point numbers)
# Examples: 1, 3.14, -5, 2.718.
x <- 10
class(x)
y <- 3.14
class(y)
a <- -5
class(a)


# 2. Character
# Character data types represent text strings.
# Examples: "hello", "R programming", "123".
name <- "John"
class(name)
my_char <- 'a'
class(my_char)


# 3. Logical
# Logical data types represent Boolean values, either TRUE or FALSE.
# Examples: TRUE, FALSE.
is_student <- TRUE
class(is_student)
is_weekday <- T
is_weekday
class(is_weekday)

# 4. Integer
# Integer data types represent whole numbers
# Examples: 1L, -5L
age <- 25L
class(age)


# 5. Complex
# Complex data types represent complex numbers with real and imaginary parts
# Example: 1 + 2i
z <- 2 + 3i
class(z)


# 6. Raw
# Raw data types represent binary data in its raw form.
# Example: as.raw(65).
raw_variable <- charToRaw("Welcome to Programm")
class(raw_variable)
char_variable <- rawToChar(raw_variable)
class(char_variable)


# 7. Factor
# Factor data types represent categorical variables with predefined levels.
# Example: factor(c("low", "medium", "high")).
gender <- factor(c("Male", "Female", "Male"))
class(gender)
example <- factor(c("low", "medium", "high"))
class(example)


# 8. Date and Time
# Date and time data types represent date and time values.
# Example: as.Date("2023-12-31"), as.POSIXct("2023-12-31 23:59:59").
birthdate <- as.Date("1990-05-15")
birthdate
class(birthdate)
timestamp <- as.POSIXct("2023-12-01 12:30:00")
timestamp
class(timestamp)


# 9. Missing Values
# Missing value data types represent missing or undefined values.
# Example: NA, NaN, NULL.


# 10. Lists
# Lists are a collection of elements, which can be of different data types.
# Example: list(1, "hello", TRUE).
my_list <- list(
  numbers = c(1, 2, 3),
  names = c("Alice", "Bob", "Charlie"),
  ages = c(25, 30, 22)
)
class(my_list)


# 11. Data Frames
# Data frames are two-dimensional structures that represent tabular data.
# Each column can be of a different data type.
# Example: data.frame(x = c(1, 2, 3), y = c("a", "b", "c")).
data <- data.frame(
  ID = c(1, 2, 3),
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 22)
)
class(data)


# 12. Vectors
# Vectors are one-dimensional arrays that can hold elements of the same data type.
# Examples: Numeric vector: c(1, 2, 3), Character vector: c("a", "b", "c").
numeric_vector <- c(1, 2, 3)
class(numeric_vector)
character_vector <- c("a", "b", "c")
class(character_vector)




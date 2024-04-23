# Functions
# 1. Built-in functions
ages <- c(22, 11, 12, 23, 24, 54, 23, 34, 56)

# min and max
min(ages)
max(ages)
range(ages)

# measure of center
mean(ages)
median(ages)
mode(ages)


# Dispersion
sd(ages)
IQR(ages)
var(ages)


# quantiles
quantile(ages, 0.25)
quantile(ages, 0.5)
quantile(ages, 0.75)

# maths
log(ages)
log10(ages)
sqrt(ages)


# 2. User defined functions / custom functions
function_name <- function (inputs) {
  # do something
}

# Basic Syntax of a Function in R
myFunction <- function(arg1, arg2, ...) {
  # Function body
  # Do something with arg1 and arg2
  # Return a value
  return(result)
}

# Simple function - add two numbers
add <- function(num1, num2) {
  total <- num1 + num2
  return(total)
}
add(20, 50)

# Function without Explicit Return
# if you don't explicitly call return(), the function will return the last expression evaluated
multiplyNumbers <- function(x, y) {
  x * y  # This will be returned automatically
}
multiplyNumbers(4, 7)

# add multiple numbers
add_mul <- function(...) {
  numbers <- c(...)
  total <- sum(numbers)
  return(total)
}

add_mul(20, 30, 50, 15)
sum(20, 30, 50)


add_new <- function(...) {
  
  # capture all argument in a list
  args <- list(...)
  
  # initialize total to zero
  total <- 0
  
  # iterate over each argument and add total
  for (num in args) {
    total <- total + num
  }
  
  # return the total
  return(total)
}

add_new(1, 3, 4, 5, 7, 8, 9)


# Function with Default Argument Values
calculateArea <- function(length = 10, width = 5) {
  area <- length * width
  return(area)
}
calculateArea()
calculateArea(length = 9)
calculateArea(width = 9)


# 3. Anonymous Functions (functions without a name)
# lambda functions
# a simple function for a short period and don't intend to reuse it elsewhere in your code
# Using sapply()
numbers <- c(1, 2, 3, 4, 5)
squared_numbers <- sapply(numbers, function(x) x^2)
squared_numbers

# Using lapply() with a List
number_list <- list(a = 1, b = 2, c = 3)
incremented_numbers <- lapply(number_list, function(x) x + 1)
incremented_numbers

# Filtering with Filter()
data <- -3:3
positive_numbers <- Filter(function(x) x > 0, data)
positive_numbers















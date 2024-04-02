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

# add two numbers
add <- function(num1, num2) {
  total <- num1 + num2
  return(total)
}
add(20, 50)

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

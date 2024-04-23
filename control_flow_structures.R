## Control Flow
# to execute code conditionally or iteratively based on certain conditions.
# if-else (if, if...else, if...else..elif)
# Loops (for, while, repeat)
# next, return, and break
library(tidyverse)


# 1. Conditional logic
num <- -9

# if statement
if(condition) {
  # code to execute if condition is TRUE
}

if (num < 0) {
  print("Negative number")
}

# if...else statement
if(condition) {
  # code to execute if condition is TRUE
} else {
  # code to execute if condition is FALSE
}

num2 <- 6 
if(num2 < 0) {
  print("Negative Number")
} else {
  print("Positive Number")
}

# ifelse() function
# Syntax
ifelse(condition, true_value, false_value)
x <- 10
ifelse(x < 0, "negative", "Postive")

# ifelse application
ages <- c(22, 11, 12, 23, 24, 54, 23, 14)
ifelse(ages < 18, "Child", "Adult")

qol_score <- c(70, 60, 43, 55, 45, 67, 80, 78, 67)
ifelse(qol_score < 50, "Poor", "Good")



## Logical Operator
# NOT
# Today is Friday (Positive Statement)
# Today is not Friday (Negative Statement) ~ Negation
2 == 2
!2 == 2
# AND
2 == 2 && 2 > 1
2 == 2 && 2 > 3
2 == 1 && 2 > 3
# OR
2 == 2 | 2 > 1
2 == 2 | 2 > 3
2 == 1 | 2 > 3


# if ... else if ... else 
# Syntax
if(condition1) {
  # code to execute if condition is TRUE
} else if(condition2) {
  # code to execute if condition1 is False and condition2 is True
} else if ( condition3) {
  # code to execute if condition1 & condition2 is False and condition3 is True
} else {
  # code to execute when all previous conditions are False
}

bmi <- 30
if(bmi < 18.5) {
  print("Underweight")
} else if(bmi >= 18.5 && bmi < 25) {
  print("Normal weight")
} else if ( bmi >= 25 && bmi < 30) {
  print("Over weight")
} else {
  print("Obese")
}

# Handling missing values with ifelse
x <- c(1, 2, NA, 3, 4, 5)
result <- ifelse(is.na(x), "Missing values", ifelse(x > 3, "Greater", "less or equeal"))
result

# if else with data frame
df <- data.frame(Name = c("John", "Micheal", "Clark", "Jack"), Age = c(24, 45, 23, 34))
df$Category <- ifelse(df$Age >= 30, "Senior", "Junior")



## 2. Loops
# for loop
for (var in seq) {
  # do something
}

# print Bangladesh five times
for (i in 1:5) {
  print("Bangladesh")
}

# print 1 to 5 
for (i in 1:5) {
  print(i)
}

# print 5 to 1
for (i in 5:1) {
  print(i)
}

# using for loop with vector
fruits <- c("Apple", "Banana", "Orange", "Jackfruit")
for (fruit in fruits) {
  print(fruit)
}

# for loop with conditions
for (i in 1:10) {
  if (i %% 2 == 0){
    print(paste(i, "is even"))
  } else{
    print(paste(i, "is odd"))
  }
}


## 3. Next statement
for (i in 1:10) {
  if(i == 3) {
    next
  }
  print(i)
}

for (i in 1:10) {
  if (i %% 2 == 0) {
    next  # Skips the rest of the loop body for even numbers
  }
  print(i)
}

# Searching for a specific condition and handling edge cases
values <- c(1, 3, 5, -1, 7, 9, 12)

for (value in values) {
  if (value == -1) {
    print("Unexpected negative value, breaking loop.")
    next  # Exit the loop if an unexpected value is found
  }
  if (value %% 2 == 0) {
    next  # Skip processing for even numbers
  }
  # Process odd positive values
  print(paste(value, "is an odd positive number."))
}


## 4. Break statement
# use to exit a loop prematurely
for (i in 1:10) {
  # Exits the loop when i equals 3
  if(i == 3) {
    break
  }
  print(i)
}

i <- 1
while(i <= 10) {
  print(i)
  if(i == 5){
    break
  }
  i <-  i + 1
}
for (i in 1:10) {
  if (i %% 2 == 0) {
    next  # Skips the rest of the loop body for even numbers
  }
  print(i)
}



## 5. While loop
while(condition) {
  # do something
}

# print 1 -10
i = 1
while (i <= 10) {
  print(i)
  i = i + 1
}

# using whilte loop with conditions
i = 1
while (i <= 10) {
  if (i %% 2 == 0) {
    print(paste(i, "is evern"))
  } else{
    print(paste(i, "is odd"))
  }
  i <- i + 1
}

# while loop with external conditions
x <- 10
while (x > 0) {
  print(x)
  x <- x - 2
}



## 6. Repeat loop
# used to repeatedly execute a block of code indefinitely until a break condition is met. 
# syntax
repeat {
  # Code to execute
  
  if (condition) {
    break  # Exit the loop if the condition is true
  }
}

i <-  1
repeat{
  print(i)
  i <-  i + 1
  if (i > 5) {
    break
  }
}

# Example: Asking for numeric input until it's provided
repeat {
  userInput <- readline(prompt="Enter a numeric value: ")
  
  # Check if the input is a numeric value
  if (is.numeric(as.numeric(userInput)) && !is.na(as.numeric(userInput))) {
    print(paste("You entered:", as.numeric(userInput)))
    break  # Exit the loop
  } else {
    print("Invalid input, please enter a numeric value.")
  }
}















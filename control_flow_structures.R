## Control Flow
# to execute code conditionally or iteratively based on certain conditions.
# if-else Statements (if, if...else, if...else..elif)
# Loops (for, while, repeat)
# Control Functions (next, return, and break)


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


# Next statement
for (i in 1:10) {
  if(i == 3) {
    next
  }
  print(i)
}


# Break statement
for (i in 1:10) {
  if(i == 3) {
    break
  }
  print(i)
}


# while loop
while(condition) {
  # do something
}

# print 1 -10
i = 1
while (i <= 10) {
  print(i)
  i = i + 1
}
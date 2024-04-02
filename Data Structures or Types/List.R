## List
# a list is a versatile data structure that can contain elements of different data types, such as vectors, matrices, data frames, or even other lists.
my_lsit <- list(
  ages <- c(22, 33, 45, 62, 34),
  smooking_statrus<- c("yes", "no", "no", "yes"),
  name <- "Ridwan"
)

# Create a list
L <- list(1, "a", TRUE, 1+3i)

# Create a list with different elements
my_list <- list(
  name = "Md Ridwan Ahmed", # character value
  age = 30, # Numeric value
  is_student = TRUE, # Logical value
  scores = c(90, 85, 92), # Numeric vector
  matrix_data =  matrix(1:6, nrow = 2, byrow = T), # Matrix
  sub_list = list(a = 10, b = c("x", "y", "z")) # Nested list
)


## Sub-setting elements of the list
# Sub-setting by name
my_list$scores 
# Sub-setting by index
my_list[[2]] 
# Sub-setting nested list element
my_list$sub_list$a 
# sub-setting a subset of list elements
my_list[2:3]


## Operations on Lists:
# length of list
length(my_list)
# Adding Elements
my_list <- c(my_list, L)


## Removing Elements:
# Remove element by index
my_list <- my_list[-3] 
# Remove all elements
my_list <- NULL             









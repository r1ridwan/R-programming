# 1. Create New Data frame
df= data.frame(col1=letters[1:5], col2=rep("r", 5),col3= 1:5)
row.names(df) <- c("abc", "abd", "acb", "cbd", "cda")


# 2. Create Data Frame
emp.data <- data.frame(
  emp_id = c (1:5), 
  emp_name = c("Rick","Dan","Michelle","Ryan","Gary"),
  salary = c(623.3,515.2,611.0,729.0,843.25), 
  
  start_date = as.Date(c("2012-01-01", "2013-09-23", "2014-11-15", "2014-05-11",
                         "2015-03-27")),
  stringsAsFactors = FALSE
)
print(emp.data)
# Get the structure of the data frame.
str(emp.data)
# Get the summary of data frame
summary(emp.data)
# Extract Specific columns
result <- data.frame(emp.data$emp_name,emp.data$salary)  
result
# Extract first two rows
emp.data[1:2,] 
# extract first 3 columns with all the rows
emp.data[,1:3] 
# Extract 3rd and 5th row with 2nd and 4th column
emp.data[c(3,5),c(2,4)] 
# New Column Addition
emp.data$dept <- c("IT","Operations","IT","HR","Finance")
emp.data.new <- emp.data
emp.data.new
# Create a Vector called DESIGNATION and Adding as a New Colum
emp.data$designation<-c("Entry level","Manager","Technical specialist","Entry level","Senior Level")
# Add new column by  cbind() function
designation <- c ("Entry level","Manager","Technical specialist","Entry level","Senior Level")
dept <- c("IT","Operations","IT","HR","Finance")
emp.table <- cbind(emp.data,dept,designation)


# Create the second data frame/ To add more rows permanently to an existing data frame, we need to bring in the new rows
# in the same structure as the existing data frame and use the rbind() function.
# In the example below we create a data frame with new rows and merge it with the existing data frame to create the final data frame.
emp.newdata <- 	data.frame(
  emp_id = c (6:8), 
  emp_name = c("Rasmi","Pranab","Tusar"),
  salary = c(578.0,722.5,632.8), 
  start_date = as.Date(c("2013-05-21","2013-07-30","2014-06-17")),
  dept = c("IT","Operations","Fianance"),
  stringsAsFactors = FALSE
)
# Bind the two data frames using rbind()
emp.finaldata <- rbind(emp.data,emp.newdata) 
print(emp.finaldata)


# 3. Create New Data frame
students_df<-data.frame(
  Subjects=c("Math", "English", "Bangla", "Science", "Sociology", "Islam" ),
  Percentage=c("90", "80", "97", "100", "87", "99")
)
# Set row names for the data frame
row.names(students_df)<-c("A", "B", "C", "D", "E", "F")
# Rename the data frame
names(students_df)<-c("Course","Score") 
# number of rows in data frame
nrow(students_df) 
# number of columns in data frame.
ncol(students_df) 
# Dimension of data frame
dim(students_df) 
# Access first row and second column of the data frame
students_df[1,2] 
# Access all the elements of the first column
students_df[,1] 


 
# Apply Function in R – apply vs lapply vs sapply vs mapply vs tapply vs rapply vs vapply
# 05. Create Data Frame
# Where the first Argument X is a data frame or matrix
# Second argument 1 indicated Processing along rows .if it is 2 then it indicated processing along the columns
# Third Argument is some aggregate function like sum, mean etc or some other user defined functions.
Age<-c(56,34,67,33,25,28)
Weight<-c(78,67,56,44,56,89)
Height<-c(165, 171,167,167,166,181)

BMI_df<-data.frame(Age,Weight,Height)
apply(BMI_df,1,sum)# row wise sum up of data frame using apply function in R
apply(BMI_df,2,sum)# column wise sum up of data frame using apply function in R
apply(BMI_df,2,mean)# column wise mean of data frame using apply function in R
#PAGE LINK: http://www.datasciencemadesimple.com/apply-function-r/

# Select elements from a data frame with the help of square brackets [ ]. By using a comma, you can indicate what to select from the rows and the columns respectively. 
# For example:
my_df[1,2] # selects the value at the first row and second column in my_df DATAFRAME.
my_df[1:3,2:4] # selects rows 1, 2, 3 and columns 2, 3, 4 in my_df DATAFRAME
my_df[1, ] #selects all elements of the first row
my_df[, 1] # Select all elements in the first column

# Sometimes you have many variables and difficult to find out the numbers of columns, in this case name of the variable is the best option
my_df[1:5, "type"] # Select first five rows and type column 


# How to Sort or Order element in a data frame? 
# In data analysis you can sort your data according to a certain variable in the dataset. In R, this is done with the help of the function order().
# order() is a function that gives you the ranked position of each element when it is applied on a variable, such as a vector for example:
a <- c(100, 10, 1000)
order(a)
a[order(a)]

  

# Basic data frame operations
# Create a data frame of boat sale data called bsale
bsale <- data.frame(name = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"),
                    color = c("black", "green", "pink", "blue", "blue", 
                              "green", "green", "yellow", "black", "black"),
                    age = c(143, 53, 356, 23, 647, 24, 532, 43, 66, 86),
                    price = c(53, 87, 54, 66, 264, 32, 532, 58, 99, 132),
                    cost = c(52, 80, 20, 100, 189, 12, 520, 68, 80, 100),
                    stringsAsFactors = FALSE)   # Don't convert strings to factors!

# Explore the bsale dataset:
head(bsale)     # Show me the first few rows
str(bsale)      # Show me the structure of the data
View(bsale)     # Open the data in a new window
names(bsale)    # What are the names of the columns?
nrow(bsale)     # How many rows are there in the data?

# Calculating statistics from column vectors
mean(bsale$age)       # What was the mean age?
table(bsale$color)    # How many boats were there of each color?
max(bsale$price)      # What was the maximum price?

# Adding new columns
bsale$id <- 1:nrow(bsale) # new colum id - value 1-10 hobe, because bsale has total 10 rows thats why nrow(bsale) function is used here
bsale$age.decades <- bsale$age / 10
bsale$profit <- bsale$price - bsale$cost

# What was the mean price of green boats?
with(bsale, mean(price[color == "green"]))

# What were the names of boats older than 100 years?
with(bsale, name[age > 100])

# What percent of black boats had a positive profit?
with(subset(bsale, color == "black"), mean(profit > 0))

# Save only the price and cost columns in a new dataframe
bsale.2 <- bsale[c("price", "cost")]

# Change the names of the columns to "p" and "c"
names(bsale.2) <- c("p", "c")

# Create a data frame called old.black.bsale containing only data from black boats older than 50 years
old.black.bsale <- subset(bsale, color == "black" & age > 50)



# CREATE A DATAFRAME OF THREE GENES WITH DIFFERENT VALUES
data = data.frame(
  C1 = c(1887.7, 9.9, 236.4), C2=c(8.1, 8.9, 8.4), C3=c(1799.6, 7.7, 220.9), 
  P1 = c(800.2, 400.1, 127.8), P2 = c(785.3, 389.3, 156.3), P3 = c(800.4, 403.9, 140.2)
)
row.names(data) <- c( "TP53", "SEZ6L", "MICALL1")

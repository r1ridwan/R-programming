## Reshape data frames between long and wide formats
# 1. Long-Wide data transformation with Melting and Casting 
# Melting function
# to transforms a data frame from a wide format to a long format
# Load packages 
install.packages("MASS")
library(MASS)
install.packages("reshape2")
library(reshape2)
install.packages("reshape")
library(reshape)

shipdata<-(head(ships,n=10))
molten.ships <- melt(shipdata, id = c("type","year")) 

# Sample data frame (wide format)
wide_data <- data.frame(
  ID = c(1, 2, 3),
  Name = c("Alice", "Bob", "Charlie"),
  Math = c(85, 78, 90),
  Science = c(75, 80, 88),
  English = c(92, 85, 78)
)
# Melting the data frame to long format
long_data <- melt(wide_data, id.vars = c("ID", "Name"), 
                  variable.name = "Subject", value.name = "Score")

# Casting function
# to transforms a data frame from a long format to a wide format
recasted.ship <- cast(molten.ships, type+year~variable,sum) 
# Sample data frame (long format)
long_data <- data.frame(
  ID = c(1, 1, 2, 2, 3, 3),
  Name = c("Alice", "Alice", "Bob", "Bob", "Charlie", "Charlie"),
  Subject = c("Math", "Science", "Math", "Science", "Math", "Science"),
  Score = c(85, 75, 78, 80, 90, 88)
)
# Casting the data frame to wide format
wide_data <- dcast(long_data, ID + Name ~ Subject, value.var = "Score")



## 2. pivot_longer() and pivot_wider() to reshape data frames
# pivot_longer() function is used to reshape data from wide to long format
# Load packages 
library(tidyr)
# Sample data frame (wide format)
wide_data <- data.frame(
  ID = c(1, 2, 3),
  Name = c("Alice", "Bob", "Charlie"),
  Math = c(85, 78, 90),
  Science = c(75, 80, 88),
  English = c(92, 85, 78)
)

# Reshaping to long format
long_data <- pivot_longer(wide_data, cols = c(Math, Science, English), 
                          names_to = "Subject", values_to = "Score")

long_data


# pivot_wider() function is used to reshape data from long to wide format
# Sample data frame (long format)
long_data <- data.frame(
  ID = c(1, 1, 2, 2, 3, 3),
  Name = c("Alice", "Alice", "Bob", "Bob", "Charlie", "Charlie"),
  Subject = c("Math", "Science", "Math", "Science", "Math", "Science"),
  Score = c(85, 75, 78, 80, 90, 88)
)

# Reshaping to wide format
wide_data <- pivot_wider(long_data, names_from = Subject, values_from = Score)
wide_data











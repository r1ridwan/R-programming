install.packages("tidyverse")
installed.packages("readr")
library(tidyverse)
library(dplyr)
library(readr)
myfile <- read_csv("U133AGNF1B.gcrma.avg.csv")

# Selecting Rows and Columns-1st Value for ROWS and 2nd Value for COLUMNS
myfile [1, 2]
myfile[, 4]
myfile[4, ]
myfile[, 4:10]
myfile[c(2, 5, 7), ]
myfile[1:3, "Adipocyte"]
myfile[120:230, "Adipocyte"]
file1<-myfile[, c("Adipocyte", "AdrenalCortex", "Adrenalgland", "Liver", "Lung")]


#EXTRACTING ROWS
myfile2<-myfile[c(1, 3, 5, 7, 22, 45), ]

#EXTRACTING ROWS WITH DIFFERENT APPROACH
myfile3<-myfile[c(1:5, 44, 200:300, 898), ]

#FILTER FUNCTION TO EXTRACT A SPECIFIC CULUMN 
myfile4<-filter(myfile, X1== "1431_at")

#FILTER MULTIPLE ROWS BY USING THEIR COLUMN NAMES
myfile5 = filter(myfile, X1 %in% c("1007_s_at", "1053_at", "117_at"))

#SELECTING ROWS USING SLICE FUNCTIONS
myfile6 <- myfile %>% slice(10:50)

#DELETING ROWS BY THEIRS ROW NUMBERS
myfile7<- myfile[-c(1:50), ]

#DEPLYR IN ACTION: Let's go through an example to see how dplyr can aid us as historians by inputting U.S. decennial census data from 1790 to 2010. 
us_state_populations_import<-read_csv("introductory_state_example.csv")
california_and_new_york_state_populations<-us_state_populations_import %>%
  filter(state %in% c("California", "New York"))
ggplot(data=california_and_new_york_state_populations, aes(x=year, y=population, color=state)) +
  geom_line() +
  geom_point()
missisipi_and_verginia_state_populations<-us_state_populations_import %>%
  filter(state %in% c("Mississippi", "Virginia"))
ggplot(data = missisipi_and_verginia_state_populations, aes(x=year, y=population, color=state))+
  geom_line() +
    geom_point()



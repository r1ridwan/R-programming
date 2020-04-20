#tidyverse istall require updated version of R---------
#install.packages("tidyverse")
library(tidyverse)
#install.packages("stringi")
library(stringi)
install.packages("stringr")
library(stringr)
#install.packages("dplyr")
library(dplyr)
#libraray and require are same function
require(gapminder)

#-------------------------------------------------------
data<- read.csv("myfile.csv")
head(data)
tail(data)
dim(data)
class(data)
data[1,3]#it will show fist row and 3rd colum data
data[100,1]#it will show 100th row and 1st colum data
data[,5]#it will show entire colum
data[11,]#it will show just entire row
data$pop #it will show the entire colum of it. 
data$continent #it will show the entire colum of it. 
#----------------------------------------------------------

data %>%
  select(continent, pop, gdpPercap)
 
temp = filter(data, country =="Afghanistan")
#if I use ! before one equal then data will be without afghanistan
temp = filter(data, country!="Afghanistan")
 #if you want to select more than one country in gapminder dataset
temp= filter(data, country %in% c("Afghanistan", "Australia"))
filter(data, year==1952)#Selecting specific year in gapminder datasets
#-------------------------
#--------------Import DAta from Excel-----------------
library(readr)
myfile <- read_csv("myfile.csv")
View(myfile)
#---------------------## DATA READING AND WRITING ##------------------------------------------
data("starwars")
#With Fixed directory
write.table(starwars, file = "new_file.csv")
head(gapminder)
#With directory
write.table(gapminder, "I:/Bioinformatics/Project 01/mygapminder_file.csv", row.names = F,
            quote = F, sep = "|" )
write_csv(gapminder, "new_file001.csv")
myfile<-read_csv(file = "new_file001.csv")

#-----------------------------STARWARS DATASETS---------------------------------------
data("starwars")
starwars %>%
  select(gender, mass, height, species) %>%
  filter(species=="Human") %>%
  na.omit() %>%
  mutate(height= height/ 100) %>% 
  mutate(BMI=mass/height^2) %>% 
  group_by(gender) %>% 
  summarise(Average_BMI=mean(BMI))
  
 head(starwars) 
tail(starwars)  
starwars[c(5,6,7,8), 5:9]  
starwars[c(6,7,8), c(2,6,9)] #starwars[c(rowss), c(colum)] 









# Loading required package
library(tidyverse)

# it will show all the existing data in R
data()

# Load Starward data
view(starwars)
data <- starwars

# Here we select 4 variable from starwars dataset. Then filter only human data.
# Then remove NA values with "na.omit()" function. Then convert height value from cm to m.
# Then create new variable for BMI with mutate functions
# Then group the data by sex and get average BMI with mutate functions
starwars %>%
  select(sex, mass, height, species) %>%
  filter(species == "Human") %>%
  na.omit() %>%
  mutate(height = height / 100) %>%
  mutate(BMI = mass / height^2) %>%
  group_by(sex) %>%
  summarise(Average = mean(BMI))


# How to use filter function

## Loading required packages and library
install.packages("dplyr")
library(dplyr)
install.packages("hflights")
library(hflights)
View(hflights)


filter(hflights, DepTime < 500 | ArrTime > 2200)->flight5

filter(hflights,Dest=="JFK" & Cancelled==1 )->flight

















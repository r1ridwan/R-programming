
#Loding requried package
library(tidyverse)

# it will show you all the existing data in R
data()

#Load Starward data
view(starwars)

# Here we select 4 variable from starwars dataset. Then filter only human data.
# Then remove NA values with "na.omit()" funciton. Then convert height value from cm to m.
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

#filter functions to filter Distance variable more than 3000
filter(hflights, Distance > 3000) -> flight1
range(flight1$Distance)

# Same function run with pipe operator, both does the same thing
hflights %>%
  filter( Distance > 2000) -> flight2


filter(hflights, UniqueCarrier %in% c("OO", "AA", "US")) -> flight1
table(flight1$UniqueCarrier)

filter(hflights, TaxiIn+TaxiOut>AirTime)->flight1
mutate(hflights,TotalTaxi=TaxiIn+TaxiOut)->flight2
filter(hflights, DepTime<500|ArrTime>2200)->flight1
filter(hflights,Dest=="JFK" & Cancelled==1 )->flight


# Using Select function to extract some specific variable and render into new dataframe 'sw'

hflights %>%
  select(Year, Month, UniqueCarrier, Distance, TaxiIn, TaxiOut) -> sw















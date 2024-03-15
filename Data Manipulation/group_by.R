library(tidyverse)
library(dplyr)
library(ggplot2)

# Load Gapminder Datasets: 
library(gapminder)

# Datasets
gapminder
iris

# group_by functions
gapminder %>%
  group_by(year) %>%
  summarise(meanlifeExp = mean(lifeExp),
            totalpop = sum(pop))

# average life expectancy and the total population in 2007 within each continent using group_by functions
gapminder %>%
  filter(year == 2007)%>%
  group_by(continent) %>%
  summarise(meanlifeExp = mean(lifeExp),
            totalpop = sum(pop))

by_continent <- gapminder %>%
  group_by(year, continent) %>%
  summarise(meanlifeExp = mean(lifeExp),
            totalpop = sum(pop))

# Find high population country in each continent 'top_n' function
gapminder %>%
  group_by(continent) %>%
  top_n(1, pop) %>%
  arrange(desc(pop))

# Find highest lifeExp for each country
gapminder_lifeExp <- gapminder %>%
  group_by(country) %>%
  top_n(1, lifeExp)



mydata2 <-iris 

# Groupby function for dataframe in R
# Mean of Sepal.Length is grouped by Species variable.
summarise_at(group_by(mydata2,Species),vars(Sepal.Length),
             funs(mean(.,na.rm=TRUE)))

# Group by function for dataframe in R using pipe operator
# Sum of Sepal.Length is grouped by Species variable with the help of pipe operator (%>%) in dplyr package. 
# As the result we will getting the sum of all the Sepal.Lengths of each species
mydata2 %>% group_by(Species) %>% 
  summarise_at(vars(Sepal.Length),
               funs(sum(.,na.rm=TRUE)))


# Groupby in R without dplyr using aggregate function:
# Sum of Sepal.Length is grouped by Species variable with the help of aggregate function in R
aggregate(mydata2$Sepal.Length, by=list(Species=mydata2$Species), FUN=sum)

# Group by function for dataframe in R using pipe operator 
# mean of Sepal.Length is grouped by Species variable with the help of pipe operator (%>%) in dplyr package. 
# As the result we will getting the mean Sepal.Length of each species
mydata2 %>% 
  group_by(Species) %>% 
  summarise_at(vars(Sepal.Length),funs(mean(.,na.rm=TRUE)))

# Groupby max in R with dplyr pipe operator %>%:
# max of Sepal.Length column is grouped by Species variable with the help of pipe operator (%>%) in dplyr package. 
# As the result we will getting the max value of Sepal.Length variable for each species
mydata2 %>% 
  group_by(Species) %>% 
  summarise_at(vars(Sepal.Length),funs(max(.,na.rm=TRUE)))

# Groupby min in R with dplyr pipe operator %>%:
# min of Sepal.Length column is grouped by Species variable with the help of pipe operator (%>%) in dplyr package. 
# As the result we will getting the min value of Sepal.Length variable for each species
mydata2 %>% 
  group_by(Species) %>% 
  summarise_at(vars(Sepal.Length),funs(min(.,na.rm=TRUE)))




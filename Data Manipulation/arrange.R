library(tidyverse)
library(dplyr)
library(ggplot2)

# Load Gapminder Datasets
library(gapminder)

# Datasets
gapminder

## Arrange function
# ascending order (lowest to highest)
gapminder %>%
  arrange(gdpPercap)
# Descending order (Highest to Lowest)
gapminder %>%
  arrange(desc(gdpPercap))
# Highest gdpPercap countries in 2007
gapminder %>%
  filter(year == 2007) %>%
  arrange(desc(gdpPercap))

# United States population 
gapminder %>%
  arrange(pop) %>%
  filter(country == "United States")

# Find Highest populations and lifeExp is greater than 50
gapminder %>%
  arrange(desc(pop)) %>%
  filter(lifeExp > 50)
gapminder %>%
  arrange(desc(pop)) %>%
  filter(lifeExp > 50, 
         country == "United States")

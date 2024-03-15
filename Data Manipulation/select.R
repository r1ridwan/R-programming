library(tidyverse)
library(dplyr)
library(ggplot2)

# Load Gapminder Datasets: 
library(gapminder)

# Datasets
gapminder

# Select function
gapminder_selected <- gapminder %>%
  select(country, year, lifeExp, pop)

# Select country, year and all the columns contain c character
gapminder %>%
  select(country, year, contains("c"))

# Select country, year and all the columns starts with l character
gapminder %>%
  select(country, year, starts_with("l"))

# Select country, year and all the columns ends with p character
gapminder %>%
  select(country, ends_with("p"))

# Remove any column from gapminder dataset
gapminder_removed <- gapminder %>%
  select(-gdpPercap)
  







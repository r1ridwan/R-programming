library(tidyverse)
library(dplyr)
library(ggplot2)

# Load Gapminder Datasets:
library(gapminder)

# Datasets
gapminder

## mutate functions
  mutate(pop = pop / 1000000)
# add new variable
gapminder %>%
  mutate(gdp = pop * gdpPercap)

# Find Highest total gdp
gapminder %>%
  mutate(gdp = pop * gdpPercap) %>%
  filter(year == 2007) %>%
  arrange(desc(gdp))










library(tidyverse)
library(dplyr)
library(ggplot2)

# Gapminder Datasets:
# it tracks economic and social indicators like life expectancy and the GDP per capital of countries over time
library(gapminder)

# Datasets
gapminder

# filter functions
gapminder %>%
  filter(year == 2007)

gapminder %>%
  filter(year == 2007 | year == 1952)

gapminder %>%
  filter(country == "United States")

gapminder %>%
  filter(year == 2007, country == "United States")

gapminder %>%
  filter(country %in% c("Bangladesh", "Afghanistan"))








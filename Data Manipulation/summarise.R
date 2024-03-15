library(tidyverse)
library(dplyr)
library(ggplot2)

# Load Gapminder Datasets: 
library(gapminder)

# Datasets
gapminder

# Summarize function
# Average life expectancy across all countries and years in the dataset
gapminder %>%
  summarise(meanlifeExp = mean(lifeExp))

# Average lifeExp in a particular year, such as 2007 and you want to find the total population in that year
gapminder %>%
  filter(year == 2007, continent == "Asia") %>%
  summarise(meanlifeExp = mean(lifeExp),
            totalpop = sum(pop))

# Median: the median represents the point in a set of numbers where half the numbers are above that point and half of the numbers are below.
gapminder %>%
  summarize (medianlifeExp = median(lifeExp))

# Max
gapminder %>%
  filter(year == 1957) %>%
  summarize(medianLifeExp = median(lifeExp),
            maxGdpPercap = max(gdpPercap))

# Lets find each continent min population, maximum lifeExp, and average GdpPercap
gapminder %>%
  group_by(continent) %>%
  summarize(min_pop = min(pop),
          max_lifeExp = max(lifeExp),
          average_gdpPercap = mean(gdpPercap))

# Install Packages
install.packages("rio")

# Load packages
library(rio)
library(tidyverse)
library(gapminder)

# Import data
data <- gapminder


## Exploring data
# 1. examine first fro rows
head(data)
head(data, 10)
head(data, n=10)

# 2. examine last few rows
tail(data)
tail(data, 10)
tail(data, n = 10)

# 3. dimension
dim(data)

# 4. number of columns
ncol(data)

# 5. number of rows
nrow(data)

# 6. check data structure 
glimpse(data)

# 7. check missing values
is.na(data)
sum(is.na(data))

# 8. check duplicated rows
duplicated(data)
sum(duplicated(data))

# 9. sampling
sample_n(data, 10)
sample_frac(data, 0.25)


## dplyr functions
# function_name(data, do something)

# 1. select
# select single column by column name
select(data, country)

# select single column by column number
select(data, 2)

# select multiple columns by column mame
select(data, country, continent, year)

# select multiple columns by column number
select(data, 1:3)
select(data, c(1, 2, 4))

# remove single column by column name
select(data, -country)

# remove single column by column number
select(data, -2)

# remove multiple columns by column name
select(data, -c(country, continent, year))

# remove multiple columns by column number
select(data, -c(1, 2, 4))

# select column by starts_with()
select(data, starts_with("c"))

# select column by ends_with()
select(data, ends_with("y"))

# select column by contains()
select(data, contains("e"))
select(data, contains("l"), contains("y"))


# 2. filter
# Equality ("==")
filter(data, country == "Bangladesh")

# Inequality("!=="
filter(data, country != "Bangladesh")

# Greater (">")
filter(data, gdpPercap > 800)

# less("<")
filter(data, gdpPercap < 800)

# greater or equal (">=")
filter(data, gdpPercap >= 800)

# less or equal ("<=")
filter(data, gdpPercap <= 800)

# Logical AND ("&")
filter(data, country == "Bangladesh" & gdpPercap >= 800)

# Logical OR ("|")
filter(data, country == "Bangladesh" | gdpPercap >= 800) 

# the "%in%" operator
filter(data, country == "India")
filter(data, country %in% c("Bangladesh", "India", "Pakistan"))


# Select and Filter
subset_data <- select(data, country, gdpPercap)
filtered_data <- filter(subset_data, country == "Bangladesh" & gdpPercap >= 800)

# chaining method ( |> or %>% ~ pipe operator)
data |>
  select(country, gdpPercap) |>
  filter(country == "Bangladesh" & gdpPercap >= 800)


# 3. arrange
# ascending order
data |>
  arrange(pop) |>
  head()

# descending order
data |>
  arrange(desc(pop)) |>
  head()


# 4. mutate
data |>
  mutate(gdp = gdpPercap * pop / 10^10) |>
  head()


# 5. rename 
data |>
  rename(population = pop) |>
  head()


# 6. group_by and summarise
data |>
  group_by(continent) |>
  summarise(mean_lifeExp = mean(lifeExp),
            median_lifeExp = median(lifeExp))








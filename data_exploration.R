# Exploring the data

# 1. check dimension
dim(gapminder)

# 2. ncol()
ncol(gapminder)

# 3. nrow()
nrow(gapminder)

# 4. colnames()
names(gapminder)

# 5. examine first few rows
head(gapminder, 20)

# examine last few rows
tail(gapminder)
tail(gapminder, n= 20)

# check data structure 
str(gapminder)
glimpse(gapminder)

# Check available data
data()
airquality
head(airquality)

# check missing values
is.na(airquality)
sum(is.na(airquality))

# check duplicated rows
duplicated(airquality)
sum(duplicated(airquality))

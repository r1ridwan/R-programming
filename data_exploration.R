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

# 6. examine last few rows
tail(gapminder)
tail(gapminder, n= 20)

# 7. check data structure 
str(gapminder)
glimpse(gapminder)

# 8. Check available data
data()
airquality
head(airquality)

# 9. check missing values
is.na(airquality)
sum(is.na(airquality))

# 10. check duplicated rows
duplicated(airquality)
sum(duplicated(airquality))

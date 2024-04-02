## Simulations - Generating Random Numbers
# runif()
?runif
runif(5)

# sample()
sample(1:10, 3)
sample(1:6, 1)

# Generating a random sample of 4 element from a vector
vector <- c('A', "B", "C", "D", "E")
sample(vector, 3)

# Normal distribution ~ rnorm()
rnorm(10, mean = 0, sd = 1)
data <- rnorm(10, mean = 0, sd = 1)
hist(data)

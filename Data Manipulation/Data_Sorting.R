# Data Sorting 

# Make up a randomly ordered vector
v = sample(101:110)
# Sort the vector
sort(v)
# Reverse sorting
sort(v, decreasing = TRUE)

# Make a data frame

df = data.frame(id = 1:4,
                weight = c(20, 27, 24,22),
                size = c("small", "large", "medium", "large"))
df 

library(plyr)
arrange(df, weight)
df[order(df$weight), ]

arrange(df, size, weight)
df[order(df$size, df$weight), ]

df[do.call(order, as.list(df)), ]


















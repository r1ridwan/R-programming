## set.seed

To generate random numbers or subsetting a portoin of sequence from millions of sequences setting
a seed is very important. In R pseudorandom number generation is done with a particular seed and if 
we know the seed we can reproduce this pseudorandom numbers. 
Example: in this example we have specified a seed for a random number generator and every time before
using this generator if we use this specific seed we will get same output, but if we don't use the seed
we will get different output
```
set.seed(123)
x = rnorm(15)
xx = rnorm(15) #it produce different output
set.seed(123)
x2 = rnorm(15) #same output
```
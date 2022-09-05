# HackerRank 10 Days of Statistics
# Central Limit Theorem 3

# A sample of 100 values
# The distribution mu = 500 and sigma = 80
# Compute the interval that covers the middle 95% of the sample mean
# In other words A, B for P(A<x<B) for z = 1.96
mu = 500
sigma = 80
z = 1.96
sample = 100

mu *= sample
sigma *= sample**0.5

a = (mu-z*sigma)/100
b = (mu+z*sigma)/100

print(round(a, 3))
print(round(b, 3))

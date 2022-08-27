# HackerRank 10 Days of Statistics
# Geometric Distribution 1

# The probability that a machine produces a defective product is 1/3
# What is the probability that the 1st defect occurs the 5th item produced?
p = 1/3
n = 5

# The formula for Geometric Distribution is:
# g(n,p) = (1-p)**(n-1)*p


def g(n, p):
    return (1-p)**(n-1)*p


print(g(n, p))

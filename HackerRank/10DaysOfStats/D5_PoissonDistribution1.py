import math
# HackerRank 10 Days of Statistics
# Poison Distribution 1

# Poisson Distribution Formula
# P(k,λ) = λ^k*e^-λ/k!
# λ: average number of successes
# k: actual number of successes


def poisson(l, k):
    return format(((l**k)*(math.e**(-l)))/math.factorial(k), '.3f')

# X follows Poisson Distribution and has a mean of 2.5
# What is the probability that the number is 5


mean = 2.5
number = 5
print(poisson(mean, number))

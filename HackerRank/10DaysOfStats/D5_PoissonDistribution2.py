import math
# HackerRank 10 Days of Statistics
# Poison Distribution 2

# Poisson Distribution Formula
# P(k,λ) = λ^k*e^-λ/k!
# λ: average number of successes
# k: actual number of successes


def poisson(l, k):
    return format(((l**k)*(math.e**(-l)))/math.factorial(k), '.3f')

# The manager of a industrial plant is planning to buy a machine
# of either type A or type B
# For each day’s operation:
# The number of repairs, X, that machine A needs is a Poisson
# random variable with mean 0.88.
# The daily cost of operating A is Ca = 160 + 40X^2
# The number of repairs, Y, that machine B needs is a Poisson
# random variable with mean 1.55.
# The daily cost of operating B is Cb = 128 + 40Y^2


a = 0.88
b = 1.55

# From theory X^2 = λ + λ^2
# Expected daily cost of A
x2 = a + a**2
print(format(160+40*x2, '.3f'))

# Expected daily cost of B
y2 = b + b**2
print(format(128+40*y2, '.3f'))

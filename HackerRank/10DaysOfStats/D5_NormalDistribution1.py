import math
import statistics
# HackerRank 10 Days of Statistics
# Normal Distribution 1

# In a certain plant, the time taken to assemble a car is a random variable, X,
# having a normal distribution with a mean of 20 hours and a standard deviation
# of 2 hours. What is the probability that a car can be assembled at this plant in:
# Less than 19.5 hours?
# Between 20 and 22 hours?

# Use of math module


def cdf(a, b, mu, sigma):
    a = 0.5*(1+math.erf((a-mu)/(sigma*math.sqrt(2))))
    b = 0.5*(1+math.erf((b-mu)/(sigma*math.sqrt(2))))
    return round(b-a, 3)


mean = 20
stdev = 2
print(cdf(0, 19.5, mean, stdev))

# 20 - 22 hours
print(cdf(20, 22, mean, stdev))

# Use of statistics module

mu = 20
sigma = 2
x = 19.5
n = statistics.NormalDist(mu, sigma)
print(round(n.cdf(x), 3))

# 20 - 22 hours

a, b = 20, 22
print(round((n.cdf(b)-n.cdf(a)), 3))

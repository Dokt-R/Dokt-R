import math
# HackerRank 10 Days of Statistics
# Normal Distribution 2

import statistics

# Students mean grade = 70, with a standard deviatio = 10
# What pecrcentage: > 80, >= 60, < 60

mu = 70
sigma = 10
n = statistics.NormalDist(mu, sigma)
print(round(100-n.cdf(80)*100, 2))
less_than_60 = round(n.cdf(60)*100, 2)
print(100-less_than_60)
print(less_than_60)

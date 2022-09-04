import statistics
# HackerRank 10 Days of Statistics
# Central Limit Theorem 1

# A large elevator can transpost 9800lbs max
# Cargo is 49 boxes
# The box weight distribution mu = 205 and sigma = 15
# What is the probability that the boxes are delivered safely?
mu = 205
sigma = 15
n = 49
max_weight = 9800
average_weight = 9800/49

mu *= n
sigma *= n**0.5

norm = statistics.NormalDist(mu, sigma)
success = norm.cdf(max_weight)
print(round(success, 5))

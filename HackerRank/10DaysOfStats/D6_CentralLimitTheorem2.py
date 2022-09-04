import statistics
# HackerRank 10 Days of Statistics
# Central Limit Theorem 2

# The number of tickets purchased follows
# The distribution mu = 2.4 and sigma = 2.0
# A few hours before game, 100 students line up to buy tickets
# If only 250 tickets are left what is the probability that all of them buy tickets.

mu = 2.4
sigma = 2.0
tickets = 250
students = 100

mu *= students
sigma *= students**0.5

norm = statistics.NormalDist(mu, sigma)
p = norm.cdf(tickets)
print(round(p, 5))

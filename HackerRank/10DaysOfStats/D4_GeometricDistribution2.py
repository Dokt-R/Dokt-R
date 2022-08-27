from math import factorial as f
# HackerRank 10 Days of Statistics
# Geometric Distribution 2

# The probability that a machine produces a defective product is 1/3
# What is the probability that the 1st defect occurs during the first 5 inspections?
p = 1/3
n = 5

# This question makes the use of Binomial distribution
# Binomial Distribution Formula:
# b(x,n,p)= n!/(x!(n-x)!)p^x*q^(n-x)


def binomial(successes, total, probability):
    return f(total)/(f(successes)*f(total-successes))*probability**successes*(1-p)**(total-successes)


atleast = 0
for x in range(1, n+1):
    atleast += binomial(x, n, p)

# Otherwise we can calculate the probability that no deffects are found
# for 5 consecutive inspections and find the answer from that. Easier solution.
no_deffect = (1-p)**5
deffect = 1 - no_deffect

print(round(atleast, 5) == round(deffect, 5))

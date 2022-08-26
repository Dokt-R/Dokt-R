from math import factorial as f
# HackerRank 10 Days of Statistics
# Binomial Distribution 2

# Binomial Distribution Formula:
# b(x,n,p)= n!/(x!(n-x)!)p^x*q^(n-x)


def binomial(successes, total, probability):
    return f(total)/(f(successes)*f(total-successes))*probability**successes*(1-p)**(total-successes)


# 12% of pistons are rejected
p = .12
q = 1-p

# From a total of 10 pistons
n = 10
pistons = 2
atleast = 0
nomore = 0
# At least 2 rejects (x=>2)
for x in range(2, n+1):
    b = binomial(x, n, p)
    atleast += b

# No more than 2 rejects (x<=2)
nomore = 0
for x in range(pistons+1):
    b = binomial(x, n, p)
    nomore += b
print(format(nomore, '.3f'))
print(format(atleast, '.3f'))

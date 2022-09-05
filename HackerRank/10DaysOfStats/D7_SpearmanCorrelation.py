# HackerRank 10 Days of Statistics
# Pearson Correlation Without Scipy

n = 10
x = list(map(float, '10 9.8 8 7.8 7.7 1.7 6 5 1.4 2'.split()))
y = list(map(float, '200 44 32 24 22 17 15 12 8 4'.split()))

rankx = sorted(x)
ranky = sorted(y)

dlist = [rankx.index(x[i])-ranky.index(y[i]) for i in range(n)]
d2 = [d**2 for d in dlist]

spearman = 1-6*sum(d2)/(n*(n**2-1))
print(round(spearman, 3))

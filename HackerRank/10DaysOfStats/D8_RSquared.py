# from sklearn.metrics import r2_score
import statistics as ss

n = 5
x = [95, 85, 80, 70, 60]
y = [85, 95, 70, 65, 70]
# print(r2_score(x, y))

# Sums
sxi = sum(x)
syi = sum(y)
sxiyi = sum([x[i]*y[i] for i in range(n)])
sxi2 = sum([x[i]**2 for i in range(n)])

# Mean
xmean = round(sxi/n, 2)
ymean = round(syi/n, 2)


sxixmean2 = sum([(xi-xmean)**2 for xi in x])
syiymean2 = sum([(yi-ymean)**2 for yi in y])
# Standart Deviation
sigmax = (sxixmean2/n)**0.5
sigmay = (syiymean2/n)**0.5

b = (n * sxiyi - sxi * syi)/(n*sxi2 - sxi**2)
a = ymean - b * xmean


def ybar(a, b, xi):
    return a + b*xi


# xi = 80
sst = syiymean2
ssr = sum([(ybar(a, b, xi)-ymean)**2 for xi in x])
sse = sum([(ybar(a, b, x[i])-y[i])**2 for i in range(n)])

r2_score = round(ssr/sst, 3)
r2 = round(1 - sse/sst, 3)

print(round(ybar(a, b, 80), 3))

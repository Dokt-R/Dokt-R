import statistics
# HackerRank 10 Days of Statistics
# Pearson Correlation

n = 10
x = list(map(float, '10 9.8 8 7.8 7.7 7 6 5 4 2'.split()))
y = list(map(float, '200 44 32 24 22 17 15 12 8 4'.split()))

# Mean
xmean = round(sum(x)/n, 2)
ymean = round(sum(y)/n, 2)

# list for variance
vx = [(xi-xmean)**2 for xi in x]
vy = [(yi-ymean)**2 for yi in y]
# Standart Deviation
sigmax = (sum(vx)/n)**0.5
sigmay = (sum(vy)/n)**0.5
# Covariance
cov = sum((x[i]-xmean)*(y[i]-ymean) for i in range(n))/n
# Pearson Correlation
pearson = cov/(sigmax*sigmay)
print(round(pearson, 3))

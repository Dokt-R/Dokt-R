# HackerRank 10 Days of Statistics
# Day 4 - Binomial Distribution 1

# First binomial exercise, solved without the use of math module
def f(n):
    if n <= 1:
        return 1
    factorial = 1
    while n > 0:
        factorial *= n
        n -= 1
    return factorial


# Boys to girls ratio is 1.09:1
b = 1.09
g = 1

# Calculate probability of a boy being born
pb = b/(b+g)
pg = 1-pb

# If there is exactly 1 child born per birth what proportion
# of families with 6 children will have at least 3 boys?
n = 6

# Calculate the Sum of probabilities because we need at least 3
# meaning that 3,4,5,or 6 kids can be boys.
sum = 0
for x in range(3, n+1):
    # Binomial Distribution Formula:
    # b(x,n,p)= n!/(x!(n-x)!)
    binomial = f(n)/(f(x)*f(n-x))*(pb**x)*(pg**(n-x))
    sum += binomial

# Print output with 3 decimals
print(format(sum, '.3f'))

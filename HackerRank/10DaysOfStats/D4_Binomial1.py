# First binomial exercise, solved without the use of math module
def f(x):
    n = 1
    while x > 0:
        n *= x
        x -= 1
    return n

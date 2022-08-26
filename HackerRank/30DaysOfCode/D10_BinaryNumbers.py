# HackerRank 30 Days of Code
# Binary Numbers

n = 6
n_binary = []
while n >= 1:
    b = n % 2
    n_binary.append(b)
    n //= 2

# If we want to print the correct binary
# n_binary.reverse()
# print("".join(map(str,n_binary)))

# Since the problem is to print the maximum consecutive 1's
count = 0
total = 0
for i in n_binary:
    if i == 1:
        count += 1
    elif i == 0:
        count = 0
    if total < count:
        total = count

print(total)

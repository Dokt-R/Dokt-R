# HackerRank 30 Days of Code
# Scope
# Compute the maximum difference between any two elements in array
class Difference:
    def __init__(self, a):
        self.__elements = a

#   The easy way
    def computeDifference(self):
        self.maximumDifference = max(self.__elements)-min(self.__elements)
#   The hard way

    def computeDifferenceFull(self):
        self.maximumDifference = 0
        compute = 0
        for value in self.__elements:
            for i in range(len(self.__elements)):
                compute = value - self.__elements[i]
                if compute > self.maximumDifference:
                    self.maximumDifference = compute


_ = 3
a = [int(e) for e in '1 2 5'.split(' ')]
d = Difference(a)
d.computeDifference()

print(d.maximumDifference)

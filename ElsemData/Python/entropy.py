'''This module is used to calculate entropies of a given time series'''

import math
import numpy as np
from scipy.signal import detrend


def shannon_info(x):
    '''This is a function to calculate Shannon Information
    (AntEntropy, i.e., opposite to Boltzmann-Gibbs-Shannon entropy)

    INPUTS:
    x = Time-Series Fragment to be analyzed

    OUTPUT:
    shannon = Shannon Information of the Time-Series Fragment

    Reference:
    Sumiyoshi Abe, "Stability of Tsallis entropy and instabilities of Renyi and normalized
    Tsallis entropies: A basis for q-exponential distributions",
    PHYSICAL REVIEW E 66, 046134 (2002)
    '''

    bins = 10
    f, xi = np.histogram(x, bins)  # Histogram Graph

    zz = np.argwhere(f)
    f = f[zz]
    f = f / sum(f)
    shannon = sum(f * np.log(f))

    return shannon


def fisher_information(x):
    '''This is a function to calculate Fisher Information

    INPUTS:
    x = Time-Series Fragment to be analyzed

    OUTPUT:
    fisher = Fisher Information of the Time-Series Fragment

    Reference:
    "Fisher information and Shannon entropy for on-line detection of transient
    signal high-values in laser Doppler flowmetry signals of healthy subjects",
    Phys. Med. Biol. 53 (2008) 5061�5076, Anne Humeau,Wojciech Trzepizur,
    David Rousseau, Francois Chapeau-Blondeau and Pierre Abraham
    '''

    bins = 10
    f, xi = np.histogram(x, bins)

    zz = np.nonzero(f)
    f = f[zz]
    f = f / sum(f)

    if len(f) == 1:
        fisher = 0
    else:
        fisher = sum((np.diff(f) ** 2) / f[0:-1])
        # If we would like to be more precise in derivetive simulation, it should be: Fisher = np.sum(((np.diff(f)/np.diff(xi))**2)/f[1:-1])
    return fisher


def tsallis_entropy(x, k=-1, q=-1):
    '''This is a function to calculate Tsallis ( non-extensive ) Entropy for specific q

    INPUTS:
    x = Time-Series Fragment to be analyzed

    q = Real-valued parameter associated to Tsallis Entropy, which quantifies
        the degree of departure from extensivity. For a physical meaning it is
        usually set to 1.8. The value q=1 is not permitted. The values q>1, q=1,
        or q<1, correspond to subextensivity, extensivity or super-extensivity,
        respectively. The parameter q behaves as a microscope for exploring
        different regions of the measure PDF : for q > 1, the more singular
        regions are amplified, while for q < 1 the less singular regions are
        accentuated.Some times a scanning of optimum q is needed,
        e.g. q = 1.1, 1.3, 1.5, 1.6, 1.8, 2, 2.5, 3

    OUTPUT:
    tsallis = Tsallis Entropy of the Time-Series Fragment

    Reference:
    Sumiyoshi Abe, "Stability of Tsallis entropy and instabilities of Renyi and
    normalized Tsallis entropies: A basis for q-exponential distributions",
    PHYSICAL REVIEW E 66, 046134 (2002)
    '''
    bins = 10
    f, xi = np.histogram(x, bins)

    # zz = find(f==0); f(zz) = eps;
    zz = np.nonzero(f)
    f = f[zz]
    f = f / sum(f)
    fq = np.power(f, q)
    tsallis = (k / (q - 1)) * (1 - sum(fq))

    return tsallis


# LOAD TIMESERIES and Get Input

#! Query Here

file = open("ts.csv")
test_signal = np.loadtxt(file, delimiter=",")
AvData = test_signal


fsize = 1024
ovl = 0
inc = round(fsize * (1 - ovl / 100))  # Increment is the "sliding" or "gliding"

# EntropyVect = [Shannon, Fisher, Tsallis], e.g.: [1,0,0] = "Calculate only Shannon Info"

# ? May Add Selection or Delete completely to save resources
EntropyVect = [1, 1, 1]

# "k" constant value, proposed: k=1, or k=1.3806503E-23 (Boltzmann constant)
#  q_in =  "q" real valued parameter, proposed: q=1.8,'

k = 1
q = 1.8

NoF = math.floor((len(AvData) - fsize) / inc) + 1

# Define Entropy Arrays for Output
ShAntEntr = np.array([])
Fi = np.array([])
FishInfo = np.array([])
TsEntr = np.array([])

for nFrame in range(NoF):  # forward frame sliding

    indx1 = (nFrame) * inc
    indx2 = indx1 + fsize
    AvDataFrame = AvData[indx1:indx2]

    # Calculate Shannon AntEntropy
    if EntropyVect[0] == 1:
        Shannon = shannon_info(AvDataFrame)
        SArray = np.repeat(Shannon, fsize)
        ShAntEntr = np.append(ShAntEntr, SArray)

    # Calculate Fisher Information
    if EntropyVect[1] == 1:
        Fisher = fisher_information(AvDataFrame)
        Fi = np.append(Fi, Fisher)
        FArray = np.repeat(Fisher, fsize)
        FishInfo = np.append(FishInfo, FArray)

    # Calculate Tsallis Entropy
    if EntropyVect[2] == 1:
        Tsallis = tsallis_entropy(AvDataFrame, k, q)
        TsArray = np.repeat(Tsallis, fsize)
        TsEntr = np.append(TsEntr, TsArray)

print(SArray, FArray, TsArray)

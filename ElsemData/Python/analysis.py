'''Implements Entropy, Hurst Exponent and Spectral Power Law Analysis on a given
time series.

REFERENCES:
Shannon:
    Sumiyoshi Abe, "Stability of Tsallis entropy and instabilities of Renyi and
    normalized Tsallis entropies: A basis for q-exponential distributions",
    PHYSICAL REVIEW E 66, 046134 (2002)

Fisher:
    "Fisher information and Shannon entropy for on-line detection of transient
    signal high-values in laser Doppler flowmetry signals of healthy subjects",
    Phys. Med. Biol. 53 (2008) 5061�5076, Anne Humeau,Wojciech Trzepizur,
    David Rousseau, Francois Chapeau-Blondeau and Pierre Abraham

Tsallis:
    Sumiyoshi Abe, "Stability of Tsallis entropy and instabilities of Renyi and
    normalized Tsallis entropies: A basis for q-exponential distributions",
    PHYSICAL REVIEW E 66, 046134 (2002)
'''

import math
import pandas as pd
import numpy as np
# from scipy.signal import detrend
# from scipy.stats import linregress, spearmanr


# ? Probably delete self.bins
# ? Argwhere and Nonzero
class SeriesAnalysis:
    '''A time series class that contains all the analysis work.'''

    def __init__(self, path, channels=None):
        '''
        k = Constant value, proposed: k=1, or k=1.3806503E-23 (Boltzmann constant)

        q = Real-valued parameter associated to Tsallis Entropy, which quantifies
        the degree of departure from extensivity. For a physical meaning it is
        usually set to 1.8. The value q=1 is not permitted. The values q>1, q=1,
        or q<1, correspond to subextensivity, extensivity or super-extensivity,
        respectively. The parameter q behaves as a microscope for exploring
        different regions of the measure PDF : for q > 1, the more singular
        regions are amplified, while for q < 1 the less singular regions are
        accentuated.Some times a scanning of optimum q is needed,
        e.g. q = 1.1, 1.3, 1.5, 1.6, 1.8, 2, 2.5, 3
        '''
        # Avoids analysing the whole document if user selects columns
        if channels is None:
            channels = list(range(2, 8))
        elif isinstance(channels, list) is False:
            channels = [channels]

        self.df = pd.read_csv(path, usecols=channels)

        # The "window" used for the analysis of the series
        self.frame_size = 1024
        self.overlap = 0
        # Increment is the "sliding" or "gliding"
        self.inc = round(self.frame_size * (1 - self.overlap / 100))

        self.k = 1
        self.q = 1.8

        self.total_frames = math.floor((len(self) - self.frame_size) /
                                       self.inc) + 1

    def __len__(self):
        '''Returns total number of rows in the dataset'''
        return len(self.df.index)

    def entropy_analysis(self):
        '''OUTPUT:
        shannon = Calculates Shannon Information (AntEntropy, i.e., opposite to
                    Boltzmann-Gibbs-Shannon entropy)
        fisher = Calculates Fisher Information
        tsallis = calculates Tsallis(non-extensive) Entropy for specific q and
                    returns Escort Tsallis Entropy
        '''
        entropies = pd.DataFrame()
        for frame_number in range(self.total_frames):  # forward frame sliding

            start = frame_number * self.inc
            stop = start + self.frame_size
            frame = self.df[start:stop]
            for column in frame:
                hist, bin_edges = np.histogram(
                    frame[column])  # Histogram Graph
                index_array = np.nonzero(hist)
                hist = hist[index_array]
                hist = hist / sum(hist)
                shannon = sum(hist * np.log(hist))
                # Fisher Entroy
                if len(hist) == 1:
                    fisher = 0
                else:
                    fisher = sum((np.diff(hist) ** 2) / hist[0:-1])
                    # If we would like to be more precise in derivetive simulation,
                    # it should be: Fisher = np.sum(((np.diff(f)/np.diff(xi))**2)/f[1:-1])

                # Tsallis Entropy
                f_powerq = np.power(hist, self.q)
                tsallis = (self.k / (self.q - 1)) * (1 - sum(f_powerq))

                results = pd.DataFrame(
                    {'channel': np.repeat(column, self.frame_size),
                     'shannon': np.repeat(-shannon, self.frame_size),
                     'fisher': np.repeat(fisher, self.frame_size),
                     'tsallis': np.repeat(tsallis, self.frame_size)},
                    index=list(i for i in range(start, stop)))
                entropies = pd.concat([entropies, results])
        # Sort each channel by index
        entropies = entropies.sort_values(by=['channel'], kind='mergesort')
        return entropies


CSV_PATH = 'C:\\Users\\Doktar\\Desktop\\git\\Dokt-R\\ElsemData\\RawData\\A2020001.csv'
analysis = SeriesAnalysis(CSV_PATH)
entropies = analysis.entropy_analysis()

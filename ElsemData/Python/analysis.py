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

from scipy.stats import linregress
import math
import pandas as pd
import numpy as np
from scipy.stats import linregress, spearmanr
# from scipy.signal import detrend


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

    def frame_values(self, frame_number):
        '''Return the values of at a given frame window
        E.g: Series[1024:2048]'''
        start = frame_number * self.inc
        stop = start + self.frame_size
        return self.df[start:stop], start, stop

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
            frame, start, stop = self.frame_values(frame_number)
            for column in frame.columns:
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
                # analysis_tuple = (['shannon', shannon], [])
                # results = pd.DataFrame(
                #     {'channel':  column,
                #      'shannon': -shannon,
                #      'fisher': fisher,
                #      'tsallis': tsallis},
                #     index=[list(i for i in range(start, stop))])
                frame_analysis = {
                    'shannon': -shannon,
                    'fisher': fisher,
                    'tsallis': tsallis}
                results = pd.DataFrame()
                for key, value in frame_analysis.items():
                    if column not in results.columns:
                        print(True)
                    # if results is None:
                    results = pd.DataFrame(
                        {'analysis': key,
                            column: value},
                        index=[start, stop])
                    # else:
                    #     next_analysis = pd.DataFrame(
                    #         {'analysis': key,
                    #          column: value},
                    #         index=[start, stop])
                    entropies = pd.concat([entropies, results])
                # entropies = pd.concat([entropies, results])
            if frame_number == 1:
                break
        # Sort each channel by index
        # entropies = entropies.sort_values(by=['channel'], kind='mergesort')
        return entropies

    def rs1(self, ts, n):
        '''
        RS1 calculates Rescale Range Analysis (by Hurst) of a time-series on a
        given vector of scales

        INPUTS:
        ts            is the input time-series vector
        n             is the vector of scales on which analysis will be performed
              is the a-priori knowledge about the model (fBm or fGn) which TS fillows

        OUTPUTS:
        log_n         is the vector of scales' logarithms
        log_rs        is the vector of mean R/S's logarithms
        '''

        N = len(ts)

        log_n = []
        log_RS = []

        for m in np.arange(0, len(n)):
            iters = np.floor(N / n[m])  # (FIX)
            if iters != 0:
                r = []
                s = []
                rs = []
                for nFrame in np.arange(0, iters):
                    indx1 = int((nFrame) * n[m])
                    indx2 = int((nFrame) * n[m] + n[m])
                    x = ts[indx1:indx2]
                    xm = np.mean(x)
                    # mean of time series

                    yN = x - xm
                    # mean-adjusted time series

                    ynN = np.cumsum(yN)  # cumulative deviate series

                    xnm = xm
                    # Ayto den to ksekatharizei stoy Eftaxia to
                    # kai einai kapws diaforetiko stoy Qian

                    r.append(max(ynN) - min(ynN))

                    s.append(math.sqrt(((sum((x - xnm) ** 2))) / n[m]))

                    if s[int(nFrame)] != 0:
                        rs.append(r[int(nFrame)] / s[int(nFrame)])
                    else:
                        s[int(nFrame)] = np.finfo(float).eps
                        rs.append(r[int(nFrame)] / s[int(nFrame)])

            log_RS.append(
                math.log2(np.mean(np.real(rs)))
            )  # Log, Log10, Log2 may be used as equivalents
            # Log, Log10, Log2 may be used as equivalents
            log_n.append(math.log2(n[m]))
        return log_RS, log_n

    def hurst_analysis(self):
        # Get maximum scale value, if desirable to have one

        Sc_max = 256
        Sc_max_pow2 = math.floor(math.log2(Sc_max))

        n_w = math.floor((len(self) - self.frame_size) / self.inc) + 1

        h = []
        log_a = []
        rr = []

        # New Code
        for frame_number in range(self.total_frames):
            frame = self.frame_values(frame_number)
            for column in frame:
                n = column.shape()  # ?
                n_max = math.floor(math.log2(n))

                n_arr = []
                for j in np.arange(1, min(Sc_max_pow2, n_max) - 1):
                    n_arr.append(2 ** (j + 2))

                # Calculation of R/S
                log_rs, log_n = rs1(column, n_arr)

                # Linear Fit of R/S

                # This is the one that gives correct estimates of r^2
                if sum(np.isnan(log_rs)) > 0:  # ! Needs Fix for NaN Values
                    # H[ind1:ind2] = NaN.*ones(w,1)
                    # log_a[ind1:ind2] = NaN.*ones(w,1)
                    # rr[ind1:ind2] = NaN.*ones(w,1)
                    print("F")
                else:
                    p = np.polyfit(log_n, log_rs, 1)
                    # p_val = np.polyval(p, log_rs)
                    # log_rs_eval = np.polyval(p, log_n)
                    r2 = linregress(log_n, log_rs).rvalue ** 2
                    while ind1 < ind2:
                        h.append(p[0])
                        log_a.append(p[1])
                        rr.append(r2)
                        ind1 += 1

        # Original code (Will be deleted)
        for q in np.arange(1, n_w + 1):
            ind1 = (q - 1) * self.inc
            ind2 = (q - 1) * self.inc + self.frame_size
            X_tot_w = self.df[ind1:ind2]
            n = len(X_tot_w)  # length of time series
            n_max = math.floor(math.log2(n))

            n_arr = []
            for j in np.arange(1, min(Sc_max_pow2, n_max) - 1):
                n_arr.append(2 ** (j + 2))

            # Calculation of R/S
            log_RS, log_n = rs1(X_tot_w, n_arr)

            # Linear Fit of R/S

            p = np.polyfit(log_n, log_RS, 1)
            # p_val = np.polyval(p, log_RS)
            # log_RS_eval = np.polyval(p, log_n)
            r2 = linregress(log_n, log_RS).rvalue ** 2
            while ind1 < ind2:
                h.append(p[0])
                log_a.append(p[1])
                rr.append(r2)
                ind1 += 1


if __name__ == '__main__':
    CSV_PATH = 'C:\\Users\\Doktar\\Desktop\\git\\Dokt-R\\ElsemData\\RawData\\A2020001.csv'
    analysis = SeriesAnalysis(CSV_PATH, ['ch1', 'ch2'])
    entropies = analysis.entropy_analysis()
    print(entropies)

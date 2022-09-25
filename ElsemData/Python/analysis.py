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


Wave Bases Functions:
    Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
    University of Colorado, Program in Atmospheric and Oceanic Sciences.
    This software may be used, copied, or redistributed as long as it is not
    sold and this copyright notice is reproduced on each copy made.  This
    routine is provided as is without any express or implied warranties
    whatsoever.
    
    Notice: Please acknowledge the use of this program in any publications:
    ``Wavelet software was provided by C. Torrence and G. Compo,
    and is available at URL: http://paos.colorado.edu/research/wavelets/''.

    Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
    Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.

    Please send a copy of such publications to either C. Torrence or G. Compo:
    Dr. Christopher Torrence               Dr. Gilbert P. Compo
    Advanced Study Program                 NOAA/CIRES Climate Diagnostics Center
    National Center for Atmos. Research    Campus Box 449
    P.O. Box 3000                          University of Colorado at Boulder
    Boulder CO 80307--3000, USA.           Boulder CO 80309-0449, USA.
    E-mail: torrence@ucar.edu              E-mail: gpc@cdc.noaa.gov
'''

import math
import pandas as pd
import numpy as np
from scipy.stats import linregress, spearmanr
# from scipy.signal import detrend
from time import perf_counter as t

# Will Delete after done with profiling
from line_profiler import LineProfiler
import cProfile
import pstats


starttime = t()
# ? Probably delete self.bins
# ? Argwhere and Nonzero


class SeriesAnalysis():
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

    def power_array(self, scale):
        '''Calculates power of 2 array based on scale and series size'''
        scale_2 = math.floor(math.log2(scale))
        n_max = math.floor(math.log2(len(self)))
        power = pd.Series([
            2**(j+2) for j in range(1, min(scale_2, n_max)-1)
        ])
        return power

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
                results = pd.DataFrame(
                    {'channel':  column,
                     'shannon': -shannon,
                     'fisher': fisher,
                     'tsallis': tsallis},
                    index=[start, stop-1])
                entropies = pd.concat([entropies, results])
        # Sort each channel by index
        entropies = entropies.sort_values(by=['channel'], kind='mergesort')
        return entropies

    def rra(self, frame, scales):
        '''
        Calculates Rescale Range Analysis (by Hurst) of a time-series on a
        given vector of scales

        INPUTS:
        frame       is the input time-series vector
        scales      is the vector of scales on which analysis will be performed

        OUTPUTS:
        log_n         is the vector of scales' logarithms
        log_rs        is the vector of mean R/S's logarithms
        '''

        # Time series model fBm. If fGn delete this line
        frame = np.diff(frame)

        log_n = []
        log_rs = []
        for scale in scales:
            iters = len(frame)//scale
            r = []
            s = []
            rs = []

            for i in range(iters):
                indx1 = i*scale
                indx2 = i*scale + scale
                subseries = frame[indx1:indx2]

                subseries_mean = subseries.sum()/subseries.size
                mean_adj_series = subseries - subseries_mean
                adjusted_squared = mean_adj_series ** 2

                # cumulative deviate series
                cum_dev = mean_adj_series.cumsum()

                r.append(max(cum_dev)-min(cum_dev))
                s.append((adjusted_squared.sum() / scale)**0.5)

                if s[i] == 0:
                    s[i] = np.finfo(float).eps
                rs.append(r[i] / s[i])

            try:
                # Log, Log10, Log2 may be used as equivalents
                log_rs.append(
                    math.log2(np.mean(np.real(rs)))
                )
                log_n.append(math.log2(scale))
            except ValueError as error:
                print(error)
                return

        return log_rs, log_n

    def hurst_analysis(self):
        '''Calculates and returns Hurst Exponent, R Squared and log(a)
        '''
        hurst_analysis = pd.DataFrame()

        # Setting maximum scale and calculating power of 2
        power_array = self.power_array(256)

        for frame_number in range(self.total_frames):
            frame, start, stop = self.frame_values(frame_number)
            for column in frame:
                # Calculation of R/S
                try:
                    s = frame[column].to_numpy()
                    log_rs, log_n = self.rra(s, power_array)

                    # Linear Fit of R/S
                    hurst, log_a = np.polyfit(log_n, log_rs, 1)
                    r_squared = linregress(log_n, log_rs).rvalue ** 2
                except TypeError as error:
                    print(f'Frame number {frame_number}: {error}')
                    # H[ind1:ind2] = NaN.*ones(w,1)
                    # log_a[ind1:ind2] = NaN.*ones(w,1)
                    # rr[ind1:ind2] = NaN.*ones(w,1)
                    hurst = log_a = r_squared = np.nan
                results = pd.DataFrame(
                    {'channel':  column,
                        'hurst': hurst,
                        'log_a': log_a,
                        'r_squared': r_squared},
                    index=[start, stop-1])
            hurst_analysis = pd.concat([hurst_analysis, results])
        # Sort each channel by index
        hurst_analysis = hurst_analysis.sort_values(
            by=['channel'], kind='mergesort')
        return hurst_analysis


if __name__ == '__main__':
    CSV_PATH = 'C:\\Users\\Doktar\\Desktop\\git\\Dokt-R\\ElsemData\\RawData\\A2020001.csv'
    analysis = SeriesAnalysis(CSV_PATH, ['ch3'])

    # f, log_a, b, rr, tt, yy, wave_j, lin_fit = powerlaw(
    #     input_signal,
    # )
    # with cProfile.Profile() as pr:
    #     a = analysis.hurst_analysis()
    # a.to_csv('hurst.csv')
    # stats = pstats.Stats(pr)
    # stats.sort_stats(pstats.SortKey.TIME)
    # stats.print_stats()


def wave_bases(mother, k, scale):
    '''WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG

    Computes the wavelet function as a function of Fourier frequency,
    used for the wavelet transform in Fourier space.
    (This program is called automatically by WAVELET)

    INPUTS:

    MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
    K = a vector, the Fourier frequencies at which to calculate the wavelet
    SCALE = a number, the wavelet scale

    OUTPUTS:

    DAUGHTER = a vector, the wavelet function
    FOURIER_FACTOR = the ratio of Fourier period to scale
    COI = a number, the cone-of-influence size at the scale
    DOFMIN = a number, degrees of freedom for each point in the wavelet power
    (either 2 for Morlet and Paul, or 1 for the DOG)
    '''

    mother = mother.upper()

    n = len(k)

    if mother == "MORLET":  # Morlet
        k0 = 6
        kTrue = np.where(k > 0.0, 1, 0)  # [1,0]  Boolean Array
        expnt = -((scale * k - k0) ** 2) / 2 * kTrue
        norm = (
            math.sqrt(scale * k[1]) * (math.pi ** (-0.25)) * math.sqrt(n)
        )  # total energy=N   [Eqn(7)]

        daughter = np.exp(expnt) * norm
        daughter = daughter * kTrue  # Heaviside step function
        fourier_factor = (4 * math.pi) / (
            k0 + math.sqrt(2 + k0 ** 2)
        )  # Scale-->Fourier [Sec.3h]
        coi = fourier_factor / math.sqrt(2)  # Cone-of-influence [Sec.3g]
        dofmin = 2  # Degrees of freedom
    elif mother == "PAUL":  # PA
        m = 4
        kTrue = np.where(k > 0.0, 1, 0)  # [1,0]  Boolean Array
        expnt = -(scale * k) * kTrue
        arr = np.arange(2, 2 * m)  # Array used in Prod
        norm = (
            math.sqrt(scale * k[1])
            * (2 ** m / math.sqrt(m * math.prod(arr)))
            * math.sqrt(n)
        )
        daughter = norm * ((scale * k) ** m) * np.exp(expnt)
        daughter = daughter * kTrue  # Heaviside step function
        fourier_factor = 4 * math.pi / (2 * m + 1)
        coi = fourier_factor * math.sqrt(2)
        dofmin = 2
    elif mother == "DOG":  # DOG
        m = 2
        expnt = -((scale * k) ** 2) / 2.0
        norm = math.sqrt(scale * k[1] / math.gamma(m + 0.5)) * math.sqrt(n)
        daughter = -norm * (1j ** m) * ((scale * k) ** m) * np.exp(expnt)
        fourier_factor = 2 * math.pi * math.sqrt(2.0 / (2 * m + 1))
        coi = fourier_factor / math.sqrt(2)
        dofmin = 1

    return daughter, fourier_factor, coi, dofmin


def wavelet2(y, dt, pad=0, dj=-1, s0=-1, j1=-1, mother=-1, sp_ty=-1):
    '''WAVELET  1D Wavelet transform with optional singificance testing

    wavelet2(Y,DT,PAD,DJ,S0,J1,mother,sp_ty)


    Custom version of wavelet.m that permits either 'log' or 'linear' scales
    spacing.

    In case of "inear" scales spacing, it is proposed to use DJ=1,
    S0=irrelevant  (e.g. = 1) and J1 = number of scales


    Computes the wavelet transform of the vector Y (length N),
    with sampling rate DT.

    By default, the Morlet wavelet (k0=6) is used.
    The wavelet basis is normalized to have total energy=1 at all scales.

    INPUTS:
    Y = the time series of length N.
    DT = amount of time between each Y value, i.e. the sampling time.

    OUTPUTS:
    WAVE is the WAVELET transform of Y. This is a complex array
    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
    The WAVELET power spectrum is ABS(WAVE)^2.
    Its units are sigma^2 (the time series variance).


    OPTIONAL INPUTS:

    *** Note *** setting any of the following to -1 will cause the default
    value to be used.

    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
    N up to the next higher power of 2. This prevents wraparound
    from the end of the time series to the beginning, and also
    speeds up the FFT's used to do the wavelet transform.
    This will not eliminate all edge effects (see COI below).

    DJ = the spacing between discrete scales. Default is 0.25.
    A smaller will give better scale resolution, but be slower to plot.

    S0 = the smallest scale of the wavelet.  Default is 2*DT.

    J1 = the of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
    to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.


    OPTIONAL OUTPUTS:

    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
    to the SCALEs.

    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
    where J1+1 is the total of scales.

    COI = if specified, then return the Cone-of-Influence, which is a vector
    of N points that contains the maximum period of useful information
    at that particular time.
    Periods greater than this are subject to edge effects.
    This can be used to plot COI lines on a contour plot by doing:
    IDL>  CONTOUR,wavelet,time,period
    IDL>  PLOTS,time,coi,NOCLIP=0
    '''
    n1 = len(y)

    # y[j * step + window], dt, pad, sc_sp, 2 *
    #         dt, nr_sc - 1, mother, sp_ty

    if s0 == -1:
        s0 = 2 * dt
    if dj == -1:
        dj = 1.0 / 4.0
    if j1 == -1:
        j1 = np.fix((math.log(n1 * dt / s0) / math.log(2)) / dj)

    # Construct time series to analyze, pad if necessary
    x = y - np.mean(y)
    if pad == 1:
        base2 = np.fix(math.log(n1) / math.log(2) +
                       0.4999)  # power of 2 nearest to N
        z = np.zeros(2 ** (int(base2) + 1) - n1)
        x = np.hstack([x, z])

    n = len(x)

    # Construct wavenumber array used in transform [Eqn(5)]
    k = np.arange(1, np.fix(n / 2) + 1)
    k = k * ((2 * math.pi) / (n * dt))
    kk = np.flip(-k, 0)
    k = np.append(k, kk[1:])
    k = np.insert(k, [0], 0)

    # Compute FFT of the (padded) time series
    f = np.fft.fft(x)  # [Eqn(3)]

    # Construct SCALE array & empty PERIOD & WAVE arrays
    if sp_ty == "log":
        scale = s0 * 2 ** (np.arange(0, j1 + 1) * dj)
    else:
        scale = np.arange(1, j1 + 1 + dj, dj) * dt
    period = scale
    # Define the Wavelet Array
    # Should be Complex type
    wave = np.zeros((int(j1 + 1), n), dtype=np.complex_)

    # Loop through all scales and compute transform
    for a1 in np.arange(0, j1 + 1):
        a1 = int(a1)
        daughter, fourier_factor, coi, dofmin = wave_bases(
            mother, k, scale[a1])
        wave[a1] = np.fft.ifft(f * daughter)  # Wavelet transform[Eqn(4)]

    period = fourier_factor * scale
    coi = (
        coi
        * dt
        * np.arange(
            1,
            ((n1 + 1) / 2 - 1),
        )
    )  # COI [Sec.3g]
    coi = np.insert(coi, 0, 1e-5)
    coi = np.append(coi, np.flip(coi))
    wave = wave[:, 0:n1]  # get rid of padding before returning
    return wave, period, scale, coi


def powerlaw(y):
    '''Power Law fitting function, using the wavelet.m function found at
    http://paos.colorado.edu/research/wavelets/ for wavelet transform.
    It is designed to be flexible in use, providing different possibilities
    of application and can be used either on GNU-Octave or on
    Matlab
    It is property of the ????? team.

    INPUT ARGUMENTS:

    y = measured time-series [vector of real numbers]
    dt = sampling period [in seconds]
    step = rolling step for the estimation of power law fitting [in samples]
    window_size = window_size for wavelet transform calculation [in samples]
    sp_ty = spacing type, 'log', or 'linear'
    nr_sc = number of scales, default = 25
    sc_sp = scale spacing, default = 0.25 for 'log' spacing, use integer for
            linear spacing
    mother = wavelet name, use any of the wave_bases.m valid wavelet names,
            like 'MORLET' or 'PAUL' or 'DOG'

    OUTPUT ARGUMENTS:

    f = frequencies corresponding to central frequencies of scales
    log_a = vector of the log10 of factor "a" of the fitted power law
            "S(f)= a * f^b" through the time-series scanning with "step"
            and "window_size"
    b = vector of the exponents "b" of the fitted power law "S(f)= a * f^b"
        through the time-series scanning with "step" and "window_size"
    rr = square (r^2) of the Spearman correlation coefficients (showing the
        quality of the linear -on log-log scale- curve fitting)
    tt = table of all time vectors corresponding to the processed
        time-series windows [hours]
    yy = table of all signal vectors corresponding to the processed
        time-series windows
    wave_j = table of the wavelet transform coefficients for all iterations
    lin_fit = linear fit parameters P for each iteration (contains log_a and
            b, but useful for applying "polyfit.m")
    '''

    nr_sc = 25
    sc_sp = 0.25
    mother = "MORLET"
    dt = 1  # Sampling period in sec
    step = 1024  # 128; # 256;
    window_size = 1024  # 1024;
    sp_ty = "log"

    pad = 0
    segment_size = len(y)
    t = np.arange(1, segment_size)
    iters = (segment_size - window_size) / step

    b = []
    log_a = []
    rr = []
    tt = []

    for j in range(int(iters)):
        window = np.arange(0, window_size)
        wave, period, scale, coi = wavelet2(
            y[j * step + window], dt, pad, sc_sp, 2 *
            dt, nr_sc - 1, mother, sp_ty
        )
        LOG_Power_Spectrum = np.log10((sum(np.abs(wave.transpose()) ** 2)))
        f = 1 / period[0:-8]  # Rn
        log_f = np.log10(1 / period[0:-8])  # Rn
        p = np.polyfit(log_f, LOG_Power_Spectrum[0:-8], 1)  # Rn

        r, pval = spearmanr(log_f, LOG_Power_Spectrum[0:-8])  # Rn
        # Profiling -> []
        b = np.append(b, p[0])
        log_a = np.append(log_a, p[1])
        rr = np.append(rr, r ** 2)
        tt = np.append(tt, t[window_size + ((j + 1) * step) - 1] * dt / 3600)

    yy = wave_j = lin_fit = False  # ! DELETE WHEN YOU GET HERE
    return f, log_a, b, rr, tt, yy, wave_j, lin_fit

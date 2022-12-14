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

import logging
import os
import mysql.connector
from collections import namedtuple
import math
import pandas as pd
import numpy as np
from scipy.stats import linregress, spearmanr
# from scipy.signal import detrend

# Logger Settings
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s -> %(message)s')
file_handler = logging.FileHandler('analysis.log', 'w')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

daily = True


class SeriesAnalysis():
    '''A time series class that contains all the analysis work.'''

    def __init__(self, path, channels=None):
        '''
        k:
            Constant, proposed: k=1, or k=1.3806503E-23 (Boltzmann constant)

        q:
            Real-valued parameter associated to Tsallis Entropy, which
            quantifies the degree of departure from extensivity. For a physical
            meaning it is usually set to 1.8. The value q=1 is not permitted.
            The values q>1, q=1, or q<1, correspond to subextensivity,
            extensivity or super-extensivity, respectively. The parameter q
            behaves as a microscope for exploring different regions of the
            measure PDF : for q > 1, the more singular regions are amplified,
            while for q < 1 the less singular regions are accentuated.
            Some times a scanning of optimum q is needed,
            e.g. q = 1.1, 1.3, 1.5, 1.6, 1.8, 2, 2.5, 3
        '''

        self.csv = pd.read_csv(path)
        self.station = self.csv.codename[1]
        self.moments = self.csv.moment
        # Avoids analysing the whole document if user selects columns
        if channels is None:
            self.df = self.csv.iloc[:, 2:8]
        elif isinstance(channels, list) is False:
            channels = [channels]
            self.df = self.csv.loc[:, channels]
        else:
            self.df = self.csv.loc[:, channels]

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
        Fram = namedtuple('Fram', ['start', 'stop', 'frame'])
        start = frame_number * self.inc
        stop = start + self.frame_size
        frame = self.df[start:stop]
        return Fram(
            start,
            stop,
            frame
        )

    def concat_column(self, dataframe, value, frame_number):
        '''
        Creates a column frame by frame and returns the result.

        Parameters:
        dataframe: The dataframe that the values will be appended to.

        value: Value that will be appended to the start and end of frame.

        frame_number: Used to calculate the start and end of frame

        Returns:
        dataframe: Dataframe with appended values to start and end of frame
                    indexes
        '''
        start, stop = self.frame_values(frame_number)[:-1]
        temp_df = pd.Series(value, index=[start, stop-1])
        dataframe = pd.concat([dataframe, temp_df])

        return dataframe

    def append_column(self, dataframe, col_name, col_df):
        '''
        Either creates index for first column or appends new columns to a
        dataframe.

        Parameters:
        dataframe: Pandas Dataframe to be returned.

        col_name: Column name that will be created.

        col_df: The dataframe that will be appended as a new column.

        Returns:
        dataframe: New dataframe with an appended column.
        '''

        try:
            dataframe[col_name] = col_df
        except ValueError:
            dataframe = pd.concat([dataframe, col_df], axis=1)
            dataframe.rename(columns={0: col_name}, inplace=True)
            if daily is True:
                nans = pd.DataFrame(index=[86016, 86399])
                dataframe = pd.concat([dataframe, nans])

        return dataframe

    def catch_datalogger_error(self, frame_number, column):
        '''
        Returns True if a given column frame is a constant value. This is
        needed since analysis on a constant value causes mathematical errors
        such as zero division.
        '''
        start = frame_number * self.inc
        stop = start + self.frame_size
        frame = self.df[start:stop]

        if np.count_nonzero(np.isnan(frame[column])) > 0:
            if np.isnan(frame[column]).all():
                logger.error("""Data entries from %s to %s are completely missing in %s.
                            Datalogger Error.""",
                             start, stop, column)
                return True
            logger.error("""Data entries from %s to %s are partially missing in %s.
                         Datalogger Error.""",
                         start, stop, column)
            return True
        if (frame[column] - np.mean(frame[column])).all() == 0:
            logger.error("""A constant value (%s) is recorded for %s at data entries %s - %s.
                         Datalogger Error.""",
                         np.mean(frame[column]), column, start, stop)
            return True
        return False

    def power_array(self, scale):
        '''Calculates power of 2 array based on scale and series size'''
        scale_2 = math.floor(math.log2(scale))
        n_max = math.floor(math.log2(len(self)))
        power = pd.Series([
            2**(j+2) for j in range(1, min(scale_2, n_max)-1)
        ])
        return power

    def entropies(self):
        '''OUTPUT:
        shannon = Calculates Shannon Information (AntEntropy, i.e., opposite to
                    Boltzmann-Gibbs-Shannon entropy)
        fisher = Calculates Fisher Information
        tsallis = calculates Tsallis(non-extensive) Entropy for specific q and
                    returns Escort Tsallis Entropy
        '''

        fisher_df = pd.DataFrame()
        tsallis_df = pd.DataFrame()
        shannon_df = pd.DataFrame()

        for column in self.df:
            fisher_col = pd.DataFrame()
            tsallis_col = pd.DataFrame()
            shannon_col = pd.DataFrame()

            for frame_number in range(self.total_frames):
                # If the frame is a constant then no analysis can be done
                if self.catch_datalogger_error(frame_number, column) is True:
                    fisher_col = self.concat_column(
                        fisher_col, np.nan, frame_number)
                    tsallis_col = self.concat_column(
                        tsallis_col, np.nan, frame_number)
                    shannon_col = self.concat_column(
                        shannon_col, np.nan, frame_number)
                    continue

                frame = self.frame_values(frame_number).frame

                hist = np.histogram(
                    frame[column], 10)[0]

                index_array = np.nonzero(hist)
                hist = hist[index_array]
                hist = hist / sum(hist)
                shannon = sum(hist * np.log(hist))/np.log(len(index_array[0]))

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

                shannon_col = self.concat_column(
                    shannon_col, -shannon, frame_number)
                fisher_col = self.concat_column(
                    fisher_col, fisher, frame_number)
                tsallis_col = self.concat_column(
                    tsallis_col, tsallis, frame_number)

            fisher_df = self.append_column(fisher_df, column, fisher_col)
            tsallis_df = self.append_column(tsallis_df, column, tsallis_col)
            shannon_df = self.append_column(shannon_df, column, shannon_col)

        return fisher_df, tsallis_df, shannon_df

    def hurst(self):
        '''
        Calculates and returns Hurst Exponent, R Squared and log(a)
        '''

        def rra(frame, scales):
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
            for scale in scales[:-2]:
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

        h_exp_df = pd.DataFrame()
        r2_df = pd.DataFrame()

        # Setting maximum scale and calculating power of 2
        power_array = self.power_array(256)

        for column in self.df:
            h_exp_col = pd.DataFrame()
            r2_col = pd.DataFrame()
            for frame_number in range(self.total_frames):
                # If the frame is a constant then no analysis can be done
                if self.catch_datalogger_error(frame_number, column) is True:
                    h_exp_col = self.concat_column(
                        h_exp_col, np.nan, frame_number)
                    r2_col = self.concat_column(r2_col, np.nan, frame_number)
                    continue

                frame = self.frame_values(frame_number).frame

                log_rs, log_n = rra(frame[column], power_array)

                # Linear Fit of R/S
                hurst = np.polyfit(log_n, log_rs, 1)
                r_squared = linregress(log_n, log_rs).rvalue ** 2

                h_exp_col = self.concat_column(
                    h_exp_col, hurst[0], frame_number)
                r2_col = self.concat_column(r2_col, r_squared, frame_number)
            h_exp_df = self.append_column(h_exp_df, column, h_exp_col)
            r2_df = self.append_column(r2_df, column, r2_col)

        return h_exp_df, r2_df

    def powerlaw(self):
        '''
        Power Law fitting function, is partly converted from a MATLAB
        file using the wavelet.m function found at
        http://paos.colorado.edu/research/wavelets/ for wavelet transform.

        Parameters:

        daily: boolean
            True:   Adds NaN values at the end to have an actual daily format
                    of analysis.
            False:  Should be set to false if analysis is done on series that
                    do not follow a daily format.

        Returns:

        power_analysis: pandas.DataFrame
            'channel': channel analyzed,
            'b_t': b(t),
            'r_squared': correlation**2

        Notes:

        Parameter daily is set to True since the usual format of the series is
        in daily files and then uploaded to a database. False is given as an
        option in case someone wants to analyze a different size of series
        E.g. a couple of hours, weekly, monthly
        '''

        def wave_bases(vector, scale):
            '''
            Computes the wavelet function as a function of Fourier frequency,
            used for the wavelet transform in Fourier space.
            (This program is called automatically by wavelet() function.)

            Parameters:

            vector: array_like (numpy.ndarray)
                Fourier frequencies at which to calculate the wavelet.
            scale: number (numpy.float64)
                The wavelet scale.

            Returns:

            daughter: array_like (numpy.ndarray)
                A vector, the wavelet function.
            fourier_factor: float
                The ratio of Fourier period to scale.
            '''

            k_0 = 6
            k_boolean = np.where(vector > 0.0, 1, 0)  # [1,0]  Boolean Array
            exponent = -((scale * vector - k_0) ** 2) / 2 * k_boolean
            norm = (
                math.sqrt(scale * vector[1]) *
                (math.pi ** (-0.25)) * math.sqrt(len(vector))
            )  # total energy=N   [Eqn(7)]

            daughter = np.exp(exponent) * norm
            daughter = daughter * k_boolean  # Heaviside step function
            fourier_factor = (4 * math.pi) / (
                k_0 + math.sqrt(2 + k_0 ** 2)
            )  # Scale-->Fourier [Sec.3h]

            return daughter, fourier_factor

        def wavelet(frame):
            '''
            Custom version of wavelet.m

            Computes the wavelet transform of the vector Y (length N),
            with sampling rate DT = 1s.

            By default, the Morlet wavelet (k_0=6) is used.
            The wavelet basis is normalized to have total energy=1 at all
            scales.

            Parameters:

            frame:
                Time series of length self.frame_size.

            Returns:

            wave: array_like (numpy.ndarray)
                The wavelet transform of frame is a complex array.
                Float(wave) gives the wavelet amplitude,
                atan(imaginary(wave),float(wave) gives the wavelet phase.
                The wavelet power spectrum is |wave|^2.
                Its units are sigma^2 (the time series variance).

            period: array_like (numpy.ndarray)
                The vector of "Fourier" periods (in time units) that
                corresponds to the scales.
            '''
            spacing = 0.25
            scales_no = 25

            # Construct time series to analyze
            timeseries = frame - np.mean(frame)
            # Construct wavenumber array used in transform [Eqn(5)]
            k = np.arange(1, np.fix(self.frame_size / 2) + 1)
            k = k * ((2 * math.pi)/self.frame_size)
            k_flip = np.flip(-k, 0)
            k = np.append(k, k_flip[1:])
            k = np.insert(k, [0], 0)

            # Compute FFT of the time series
            fast_fourier = np.fft.fft(timeseries)  # [Eqn(3)]

            # Construct scale array & empty period & wave arrays
            scales = 2 * 2**(np.arange(0, scales_no) * spacing)

            # Define the Wavelet Array. Should be Complex type
            wave = np.zeros((scales_no, self.frame_size), dtype=np.complex_)

            # Loop through all scales and compute transform
            for i, scale in enumerate(scales):
                daughter, fourier_factor = wave_bases(k, scale)
                # Wavelet transform[Eqn(4)]
                wave[i] = np.fft.ifft(fast_fourier * daughter)

            period = fourier_factor * scales
            return wave, period

        def calculations(wave, period):
            '''
            A help function that calculates b(t) and r^2 given wave and period.

            Parameters:

            wave: array_like (numpy.ndarray)
                The wavelet transform of frame is a complex array.
                Float(wave) gives the wavelet amplitude,
                atan(imaginary(wave),float(wave) gives the wavelet phase.
                The wavelet power spectrum is |wave|^2.
                Its units are sigma^2 (the time series variance).

            period: array_like (numpy.ndarray)
                The vector of "Fourier" periods (in time units) that
                corresponds to the scales.

            Returns:

            b(t), r**2
            '''
            log_power_spectrum = np.log10(
                (sum(np.abs(wave.transpose()) ** 2)))

            log_f = np.log10(1 / period[0:-8])

            # p[0]=b(t), p[1]=log_a
            p = np.polyfit(log_f, log_power_spectrum[0:-8], 1)

            correlation = spearmanr(
                log_f, log_power_spectrum[0:-8]).correlation

            return p[0], correlation**2

        b_df = pd.DataFrame()
        r_df = pd.DataFrame()

        for column in self.df:
            b_col = pd.DataFrame()
            r_col = pd.DataFrame()
            for frame_number in range(self.total_frames):
                # If the frame is a constant then no analysis can be done
                if self.catch_datalogger_error(frame_number, column) is True:
                    b_col = self.concat_column(b_col, np.nan, frame_number)
                    r_col = self.concat_column(r_col, np.nan, frame_number)
                    continue
                frame = self.frame_values(frame_number).frame
                wave, period = wavelet(frame[column])
                b_t, r_squared = calculations(wave, period)

                b_col = self.concat_column(b_col, b_t, frame_number)
                r_col = self.concat_column(r_col, r_squared, frame_number)
            b_df = self.append_column(b_df, column, b_col)
            r_df = self.append_column(r_df, column, r_col)

        return b_df, r_df


if __name__ == '__main__':
    CSV_PATH = 'C:\\Users\\Doktar\\Desktop\\git\\Dokt-R\\ElsemData\\RawData\\A201.csv'
    analysis = SeriesAnalysis(CSV_PATH, ['ch6', 'ch2'])

    fisher, tsallis, shannon = analysis.entropies()
    hurst_exponent, hurst_r2 = analysis.hurst()
    b_t, power_r2 = analysis.powerlaw()

    print(b_t)
    print(power_r2)

    fisher.to_csv('fisher')
    tsallis.to_csv('tsallis')
    shannon.to_csv('shannon')
    hurst_exponent.to_csv('hurst_exponent')
    hurst_r2.to_csv('hurst_r2')
    b_t.to_csv('b_t')
    power_r2.to_csv('power_r2')

    # print(pd.read_csv(CSV_PATH, index_col='moment'))

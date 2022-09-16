# Main file to apply spectral power law
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


#       WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG
#
#      Computes the wavelet function as a function of Fourier frequency,
#      used for the wavelet transform in Fourier space.
#      (This program is called automatically by WAVELET)
#
#    INPUTS:
#
#       MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
#       K = a vector, the Fourier frequencies at which to calculate the wavelet
#       SCALE = a number, the wavelet scale
#       PARAM = the nondimensional parameter for the wavelet function
#
#    OUTPUTS:
#
#       DAUGHTER = a vector, the wavelet function
#       FOURIER_FACTOR = the ratio of Fourier period to scale
#       COI = a number, the cone-of-influence size at the scale
#       DOFMIN = a number, degrees of freedom for each point in the wavelet power
#                (either 2 for Morlet and Paul, or 1 for the DOG)
#
#   ----------------------------------------------------------------------------
#      Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
#      University of Colorado, Program in Atmospheric and Oceanic Sciences.
#      This software may be used, copied, or redistributed as long as it is not
#      sold and this copyright notice is reproduced on each copy made.  This
#      routine is provided as is without any express or implied warranties
#      whatsoever.
#   ----------------------------------------------------------------------------


def wave_bases(mother, k, scale, param):

    mother = mother.upper()

    n = len(k)

    if mother == "MORLET":  # Morlet
        if param == -1:
            param = 6
        k0 = param
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
    elif mother == "PAUL":  # Paul
        if param == -1:
            param = 4
        m = param
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
        if param == -1:
            param = 2
        m = param
        expnt = -((scale * k) ** 2) / 2.0
        norm = math.sqrt(scale * k[1] / math.gamma(m + 0.5)) * math.sqrt(n)
        daughter = -norm * (1j ** m) * ((scale * k) ** m) * np.exp(expnt)
        fourier_factor = 2 * math.pi * math.sqrt(2.0 / (2 * m + 1))
        coi = fourier_factor / math.sqrt(2)
        dofmin = 1
    else:
        #! This is error handling and must be fixed to suit python
        print("Mother must be one of MORLET,PAUL,DOG")
    return daughter, fourier_factor, coi, dofmin


#   WAVELET  1D Wavelet transform with optional singificance testing
#
#   wavelet2(Y,DT,PAD,DJ,S0,J1,mother,sp_ty,param)
#
#
#   Custom version of wavelet.m that permits either 'log' or 'linear' scales
#   spacing.
#
#   In case of "inear" scales spacing, it is proposed to use DJ=1,
#   S0=irrelevant #  (e.g. = 1) and J1 = number of scales
#
#
#   Computes the wavelet transform of the vector Y (length N),
#   with sampling rate DT.
#
#   By default, the Morlet wavelet (k0=6) is used.
#   The wavelet basis is normalized to have total energy=1 at all scales.
#
#
# INPUTS:
#
#    Y = the time series of length N.
#    DT = amount of time between each Y value, i.e. the sampling time.
#
# OUTPUTS:
#
#    WAVE is the WAVELET transform of Y. This is a complex array
#    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
#    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
#    The WAVELET power spectrum is ABS(WAVE)^2.
#    Its units are sigma^2 (the time series variance).
#
#
# OPTIONAL INPUTS:
#
# *** Note *** setting any of the following to -1 will cause the default
#               value to be used.
#
#    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
#         N up to the next higher power of 2. This prevents wraparound
#         from the end of the time series to the beginning, and also
#         speeds up the FFT's used to do the wavelet transform.
#         This will not eliminate all edge effects (see COI below).
#
#    DJ = the spacing between discrete scales. Default is 0.25.
#         A smaller # will give better scale resolution, but be slower to plot.
#
#    S0 = the smallest scale of the wavelet.  Default is 2*DT.
#
#    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
#        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
#
#
# OPTIONAL OUTPUTS:
#
#    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
#           to the SCALEs.
#
#    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
#            where J1+1 is the total # of scales.
#
#    COI = if specified, then return the Cone-of-Influence, which is a vector
#        of N points that contains the maximum period of useful information
#        at that particular time.
#        Periods greater than this are subject to edge effects.
#        This can be used to plot COI lines on a contour plot by doing:
#            IDL>  CONTOUR,wavelet,time,period
#            IDL>  PLOTS,time,coi,NOCLIP=0
#
#
# ----------------------------------------------------------------------------
#   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
#   University of Colorado, Program in Atmospheric and Oceanic Sciences.
#   This software may be used, copied, or redistributed as long as it is not
#   sold and this copyright notice is reproduced on each copy made.  This
#   routine is provided as is without any express or implied warranties
#   whatsoever.
#
#   Notice: Please acknowledge the use of this program in any publications:
#   ``Wavelet software was provided by C. Torrence and G. Compo,
#     and is available at URL: http://paos.colorado.edu/research/wavelets/''.
#
#   Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
#            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
#
#   Please send a copy of such publications to either C. Torrence or G. Compo:
#   Dr. Christopher Torrence               Dr. Gilbert P. Compo
#   Advanced Study Program                 NOAA/CIRES Climate Diagnostics Center
#   National Center for Atmos. Research    Campus Box 449
#   P.O. Box 3000                          University of Colorado at Boulder
#   Boulder CO 80307--3000, USA.           Boulder CO 80309-0449, USA.
#   E-mail: torrence@ucar.edu              E-mail: gpc@cdc.noaa.gov
# ----------------------------------------------------------------------------


def wavelet2(y, dt, pad=0, dj=-1, s0=-1, j1=-1, mother=-1, sp_ty=-1, param=-1):
    n1 = len(y)

    if s0 == -1:
        s0 = 2 * dt
    if dj == -1:
        dj = 1.0 / 4.0
    if j1 == -1:
        j1 = np.fix((math.log(n1 * dt / s0) / math.log(2)) / dj)
    if mother == -1:
        mother = "MORLET"
        sp_ty = "log"

    # Construct time series to analyze, pad if necessary
    x = y - np.mean(y)
    if pad == 1:
        base2 = np.fix(math.log(n1) / math.log(2) + 0.4999)  # power of 2 nearest to N
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
    wave = np.zeros((int(j1 + 1), n), dtype=np.complex_)  # Should be Complex type

    # Loop through all scales and compute transform
    for a1 in np.arange(0, j1 + 1):
        a1 = int(a1)
        daughter, fourier_factor, coi, dofmin = wave_bases(mother, k, scale[a1], param)
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


# Power Law fitting function, using the wavelet.m function found at
# http://paos.colorado.edu/research/wavelets/ for wavelet transform.
# It is designed to be flexible in use, providing different possibilities
# of application and can be used either on GNU-Octave or on
# Matlab
# It is property of the ????? team.
#
# INPUT ARGUMENTS:
#
# y = measured time-series [vector of real numbers]
# dt = sampling period [in seconds]
# step = rolling step for the estimation of power law fitting [in samples]
# window_size = window_size for wavelet transform calculation [in samples]
# sp_ty = spacing type, 'log', or 'linear'
# nr_sc = number of scales, default = 25
# sc_sp = scale spacing, default = 0.25 for 'log' spacing, use integer for
#         linear spacing
# mother = wavelet name, use any of the wave_bases.m valid wavelet names,
#          like 'MORLET' or 'PAUL' or 'DOG'
# plot_flag = flag to check plotting inside the function.
#            If = 0 do not plot
#            If >1 plot (plot_flag) number of figures (if available)
# print_shift = how many iterations to shift before checking to print #
# octave_flag = flag to check Octave run. If equals 1 then run is in Octave
#                                        If equals 0 then run is in MATLAB
# fast_flag = flag to indicate fast version of the routine, where memory
#             consuming computations (wave_j, lin_fit and yy) are ommited
#
# OUTPUT ARGUMENTS:
#
# f = frequencies corresponding to central frequencies of scales
# log_a = vector of the log10 of factor "a" of the fitted power law
#          "S(f)= a * f^b" through the time-series scanning with "step"
#          and "window_size"
# b = vector of the exponents "b" of the fitted power law "S(f)= a * f^b"
#      through the time-series scanning with "step" and "window_size"
# rr = square (r^2) of the Spearman correlation coefficients (showing the
#      quality of the linear -on log-log scale- curve fitting)
# tt = table of all time vectors corresponding to the processed
#      time-series windows [hours]
# yy = table of all signal vectors corresponding to the processed
#      time-series windows
# wave_j = table of the wavelet transform coefficients for all iterations
# lin_fit = linear fit parameters P for each iteration (contains log_a and
#           b, but useful for applying "polyfit.m")


def powerlaw(
    y,
    dt,
    step=-1,
    window_size=-1,
    sp_ty="log",
    nr_sc=-1,
    sc_sp=-1,
    mother=-1,
    r_critical=-1,
    plot_flag=-1,
    print_shift=-1,
    octave_flag=-1,
    fast_flag=-1,
):

    if fast_flag == -1:
        fast_flag = 1
    if octave_flag == -1:
        octave_flag = 1
    if print_shift == -1:
        print_shift = 50
    if plot_flag == -1:
        plot_flag = 0
    if r_critical == -1:
        r_critical = 0.95
    if mother == -1:
        mother = "MORLET"
    if sc_sp == -1:
        sc_sp = 0.25
    if nr_sc == -1:
        nr_sc = 25
    if window_size == -1:
        window_size = 128
    if step == -1:
        step = 1

    pad = 0
    segment_size = len(y)
    t = np.arange(1, segment_size)
    iters = (segment_size - window_size) / step
    rng = 60

    N_of_passes = 0
    fs = 1 / dt
    if octave_flag == 0:
        wvname = "morl"

    rng = 50
    bw = 10

    step = 1024  #! DELETE AFTER TEST

    b = []
    log_a = []
    rr = []
    tt = []

    for j in range(int(iters)):
        window = np.arange(0, window_size)
        wave, period, scale, coi = wavelet2(
            y[j * step + window], dt, pad, sc_sp, 2 * dt, nr_sc - 1, mother, sp_ty
        )
        LOG_Power_Spectrum = np.log10((sum(np.abs(wave.transpose()) ** 2)))
        f = 1 / period[0:-8]  # Rn
        log_f = np.log10(1 / period[0:-8])  # Rn
        p = np.polyfit(log_f, LOG_Power_Spectrum[0:-8], 1)  # Rn

        r, pval = spearmanr(log_f, LOG_Power_Spectrum[0:-8])  # Rn
        b = np.append(b, p[0])
        log_a = np.append(log_a, p[1])
        rr = np.append(rr, r ** 2)
        tt = np.append(tt, t[window_size + ((j + 1) * step) - 1] * dt / 3600)

    yy = wave_j = lin_fit = False  #! DELETE WHEN YOU GET HERE
    return f, log_a, b, rr, tt, yy, wave_j, lin_fit


#! Test Power Law Function

# function callPowerlaw(Input_Signal,dt,step,window_size,octave_flag)
#
# This function is based on the Power Law fitting function powerlaw
# It is designed to be flexible in use, providing different possibilities
# of application and can be used either on GNU-Octave or on
# Matlab
# It is property of the ????? team.
#
# INPUT ARGUMENTS:
#
# Input_Signal = measured time-series [vector of real numbers]
# dt = sampling period [in seconds]
# step = rolling step for the estimation of power law fitting [in samples]
# window_size = window_size for wavelet transform calculation [in samples]
# octave_flag = flag to check Octave run. If equals 1 then run is in Octave
#                                        If equals 0 then run is in MATLAB


def callPowerlaw(Input_Signal, dt, step, window_size, octave_flag):
    Ni = len(Input_Signal)  # extract # of data-series points
    sp_ty = "log"
    nr_sc = 25
    sc_sp = 0.25
    mother = "MORLET"
    plot_flag = 0
    r_critical = 0.85
    print_shift = 50
    fast_flag = 1

    f, log_a, b, rr, tt, yy, wave_j, lin_fit = powerlaw(
        Input_Signal,
        dt,
        step,
        window_size,
        sp_ty,
        nr_sc,
        sc_sp,
        mother,
        r_critical,
        plot_flag,
        print_shift,
        octave_flag,
        fast_flag,
    )

    return f, log_a, b, rr, tt, yy, wave_j, lin_fit


#! Load File and make it into an ndarray

#! Query Here

file = open("ts.csv")
input = np.loadtxt(file, delimiter=",")

dt = 1  # Sampling period in sec
step = 1024  # 128; # 256;
window_size = 1024  # 1024;
octave_flag = 0

#! Execute Code

f, log_a, b, rr, tt, yy, wave_j, lin_fit = callPowerlaw(
    input, dt, step, window_size, octave_flag
)

print(f, log_a, b, rr, tt, yy, wave_j, lin_fit)
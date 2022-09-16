import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def rs1(ts, n, ts_model):
    # RS1 calculates Rescale Range Analysis (by Hurst) of a time-series on a given vector of scales
    #
    # INPUTS:
    # ts            is the input time-series vector
    # n             is the vector of scales on which analysis will be performed
    # ts_model      is the a-priori knowledge about the model (fBm or fGn) which TS fillows
    #
    # OUTPUTS:
    # log_n         is the vector of scales' logarithms
    # log_rs        is the vector of mean R/S's logarithms

    if ts_model == 1:
        ts = np.diff(ts)

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
        log_n.append(math.log2(n[m]))  # Log, Log10, Log2 may be used as equivalents
    return log_RS, log_n


# LOAD TIMESERIES and Get Input

#! Query Here

# Test Signal
file = open("ts.csv")
test_signal = np.loadtxt(file, delimiter=",")
X_tot = test_signal

# Desirable window in # of points (e.g. "1024")

w = 1024

# Get the percentage of overlapping The overlapping leads to "sliding" or
# "gliding" the frame window along the time-series,
# if there is no overlap (ovl=0), then we have "lumping"
# Take into account overlapping to determine step increment in samples

ovl = 0  # 0 <= ovl <= 100

# Increment is the "sliding" or "gliding" step

inc = round(w * (1 - (ovl / 100)))

# Get maximum scale value, if desirable to have one

Sc_max = 256

Sc_max_pow2 = math.floor(math.log2(Sc_max))

# Model of the TS: 1 = fBm, 2 = fGn

ts_model = 1

# Just one area will be analyzed (MAY DELETE AFTER PLOT)
ar_num = 1
num = ar_num * 2

N_tot = len(X_tot)

n_w = math.floor((N_tot - w) / inc) + 1

h = []
log_a = []
rr = []

for q in np.arange(1, n_w + 1):
    ind1 = (q - 1) * inc
    ind2 = (q - 1) * inc + w
    X_tot_w = X_tot[ind1:ind2]
    n = len(X_tot_w)  # length of time series
    n_max = math.floor(math.log2(n))

    n_arr = []
    for j in np.arange(1, min(Sc_max_pow2, n_max) - 1):
        n_arr.append(2 ** (j + 2))

    ## Calculation of R/S
    log_RS, log_n = rs1(X_tot_w, n_arr, ts_model)

    ## Linear Fit of R/S

    # This is the one that gives correct estimates of r^2
    if sum(np.isnan(log_RS)) > 0:  #! Need Fix
        # H[ind1:ind2] = NaN.*ones(w,1)
        # log_a[ind1:ind2] = NaN.*ones(w,1)
        # rr[ind1:ind2] = NaN.*ones(w,1)
        print("F")
    else:
        p = np.polyfit(log_n, log_RS, 1)
        # p_val = np.polyval(p, log_RS)
        # log_RS_eval = np.polyval(p, log_n)
        r2 = linregress(log_n, log_RS).rvalue ** 2
        while ind1 < ind2:
            h.append(p[0])
            log_a.append(p[1])
            rr.append(r2)
            ind1 += 1
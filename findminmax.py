import numpy as np
import matplotlib.pyplot as plt

def from_timeseries(x,bins=50,npoly=17):
    n,bins = np.histogram(x,bins=bins)
    bin_avg = 0.5*(bins[1:] + bins[:-1])
    Fdata = -np.log(n)
    Fdata -= min(Fdata)

    xinterp = np.linspace(bin_avg.min(),bin_avg.max(),500)
    F = np.poly1d(np.polyfit(bin_avg,Fdata,npoly))
    dFdx = np.polyder(F,m=1)
    d2Fdx2 = np.polyder(F,m=2)

    # Find all minima and maxima
    A = dFdx(xinterp)
    minidx = np.where([(A[i] < 0) & (A[i + 1] > 0) for i in range(len(A) - 1) ])[0]
    maxidx = np.where([(A[i] > 0) & (A[i + 1] < 0) for i in range(len(A) - 1) ])[0]
    minxval = xinterp[minidx]
    maxxval = xinterp[maxidx]
    omegamin = d2Fdx2(minxval)
    omegamax = d2Fdx2(maxxval)

    return xinterp, F, minxval, maxxval, omegamin, omegamax

def from_profile(bins,Fdata,npoly=17):
    
    bin_avg = 0.5*(bins[1:] + bins[:-1])
    Fdata -= min(Fdata)

    xinterp = np.linspace(bin_avg.min(),bin_avg.max(),500)
    F = np.poly1d(np.polyfit(bin_avg,Fdata,npoly))
    dFdx = np.polyder(F,m=1)
    d2Fdx2 = np.polyder(F,m=2)

    # Find all minima and maxima
    A = dFdx(xinterp)
    minidx = np.where([(A[i] < 0) & (A[i + 1] > 0) for i in range(len(A) - 1) ])[0]
    maxidx = np.where([(A[i] > 0) & (A[i + 1] < 0) for i in range(len(A) - 1) ])[0]
    minxval = xinterp[minidx]
    maxxval = xinterp[maxidx]
    omegamin = d2Fdx2(minxval)
    omegamax = d2Fdx2(maxxval)

    return xinterp, F, minxval, maxxval, omegamin, omegamax

def state_bounds(x,F,minxval,maxxval):

    minidx = np.array([ np.where(x == minxval[i])[0][0] for i in range(minxval.shape[0]) ])
    maxidx = np.array([ np.where(x == maxxval[i])[0][0] for i in range(maxxval.shape[0]) ])

    min_state_bounds = []
    for i in range(minxval.shape[0]):
        left_min_bound = x.min()
        right_min_bound = x.max()
        for j in range(x[:minidx[i]].shape[0]):
            deltaF = F_fit(x[minidx[i] - j]) - F_fit(x[minidx[i]])
            if deltaF >= 0.3:
                left_min_bound = x[minidx[i] - j]
                break
        for j in range(x[minidx[i]:].shape[0]):
            deltaF = F_fit(x[minidx[i] + j]) - F_fit(x[minidx[i]])
            if deltaF >= 0.3:
                right_min_bound = x[minidx[i] + j]
                break
        min_state_bounds.append([left_min_bound,right_min_bound])

    max_state_bounds = []
    for i in range(minxval.shape[0]):
        left_min_bound = x.min()
        right_min_bound = x.max()
        for j in range(x[:maxidx[i]].shape[0]):
            deltaF = abs(F_fit(x[maxidx[i] - j]) - F_fit(x[maxidx[i]]))
            if deltaF >= 0.3:
                left_min_bound = x[maxidx[i] - j]
                break
        for j in range(x[maxidx[i]:].shape[0]):
            deltaF = abs(F_fit(x[maxidx[i] + j]) - F_fit(x[maxidx[i]]))
            if deltaF >= 0.3:
                right_min_bound = x[maxidx[i] + j]
                break
        max_state_bounds.append([left_min_bound,right_min_bound])

    return min_state_bounds, max_state_bounds



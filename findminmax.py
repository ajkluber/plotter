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

    return xinterp, F, minxval, maxxval

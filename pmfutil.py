import numpy as np

def pmf1D(xvst,bins=50):
    """Histogram timeseries to get 1D pmf"""
    n,bins = np.histogram(xvst,bins=bins)
    mid_bin = 0.5*(bins[1:] + bins[:-1])
    Fdata = -np.log(n)
    Fdata -= min(Fdata)
    return mid_bin,Fdata

def interpolate_profile(mid_bin,Fdata,npoly=20):
    """Interpolate 1D profile with polynomial"""
    xinterp = np.linspace(mid_bin.min(),mid_bin.max(),500)
    F = np.poly1d(np.polyfit(mid_bin,Fdata,npoly))
    return xinterp, F

def extrema_from_profile(xinterp,F):
    """Find extrema of a 1D free energy profile"""
    dFdx = np.polyder(F,m=1)
    A = dFdx(xinterp)
    minidx = np.where([(A[i] < 0) & (A[i + 1] > 0) for i in range(len(A) - 1) ])[0]
    maxidx = np.where([(A[i] > 0) & (A[i + 1] < 0) for i in range(len(A) - 1) ])[0]
    return minidx,maxidx

def extrema_from_timeseries(xvst,bins=50,npoly=20):
    """Find the minima(maxima) along 1D free energy profile from timeseries"""
    mid_bin, Fdata = pmf1D(xvst,bins=bins)
    xinterp, F = interpolate_profile(mid_bin,Fdata,npoly=20)
    minidx, maxidx = extrema_from_profile(xinterp,F)
    return minidx, maxidx

def second_deriv_from_profile(xinterp,F):
    """Calculate second derivative at extrema"""
    dFdx = np.polyder(F,m=1)
    d2Fdx2 = np.polyder(F,m=2)
    minidx, maxidx = extrema_from_profile(xinterp,F)
    omegamin = d2Fdx2(xinterp[minidx])
    omegamax = d2Fdx2(xinterp[maxidx])
    return omegamin, omegamax

def state_bounds_from_profile(xinterp,F,threshold=0.3):
    """Find boundaries of each extrema state along 1D profile"""

    minidx, maxidx = extrema_from_profile(xinterp,F)

    min_state_bounds = []
    for i in range(minxval.shape[0]):
        left_min_bound = xinterp.min()
        right_min_bound = xinterp.max()
        for j in range(xinterp[:minidx[i]].shape[0]):
            deltaF = F_fit(xinterp[minidx[i] - j]) - F_fit(xinterp[minidx[i]])
            if deltaF >= threshold:
                left_min_bound = x[minidx[i] - j]
                break
        for j in range(xinterp[minidx[i]:].shape[0]):
            deltaF = F_fit(xinterp[minidx[i] + j]) - F_fit(xinterp[minidx[i]])
            if deltaF >= threshold:
                right_min_bound = xinterp[minidx[i] + j]
                break
        min_state_bounds.append([left_min_bound,right_min_bound])

    max_state_bounds = []
    for i in range(minxval.shape[0]):
        left_min_bound = xinterp.min()
        right_min_bound = xinterp.max()
        for j in range(xinterp[:maxidx[i]].shape[0]):
            deltaF = abs(F_fit(xinterp[maxidx[i] - j]) - F_fit(xinterp[maxidx[i]]))
            if deltaF >= threshold:
                left_min_bound = xinterp[maxidx[i] - j]
                break
        for j in range(xinterp[maxidx[i]:].shape[0]):
            deltaF = abs(F_fit(xinterp[maxidx[i] + j]) - F_fit(xinterp[maxidx[i]]))
            if deltaF >= threshold:
                right_min_bound = xinterp[maxidx[i] + j]
                break
        max_state_bounds.append([left_min_bound,right_min_bound])

    return min_state_bounds, max_state_bounds

def assign_state_labels(min_bounds,max_bounds):
    """Label extrema along profile"""
    min_labels = []
    leftbounds = [ min_bounds[i][0] for i in range(len(min_bounds)) ]
    mina = min(leftbounds)
    maxa = max(leftbounds)
    for i in range(len(min_bounds)): 
        a = min_bounds[i][0]
        b = min_bounds[i][1]
        if a == mina:
            state = "U"
        elif a == maxa:
            state = "N"
        else:
            state = "I%d" % (counter)
            counter += 1
        min_labels.append(state)
    max_labels = []
    counter = 1
    for i in range(len(max_bounds)): 
        a = max_bounds[i][0]
        b = max_bounds[i][1]
        state = "TS%d" % (counter)
        counter += 1
        max_labels.append(state)
    return min_labels, max_labels

def save_state_bounds(coord_name,min_bounds,max_bounds,min_labels,max_labels): 
    """Write state bounds to file. Label intermediates"""
    with open("%s_state_bounds.txt" % coord_name,"w") as fout:
        for i in range(len(min_bounds)): 
            state_string = "%s  %e  %e\n" % (min_labels[i],min_bounds[i][0],min_bounds[i][1])
            fout.write(state_string)

        for i in range(len(max_bounds)): 
            state_string = "%s  %e  %e\n" % (max_labels[i],max_bounds[i][0],max_bounds[i][1])
            fout.write(state_string)


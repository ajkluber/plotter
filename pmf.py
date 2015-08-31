import argparse
import numpy as np
import matplotlib.pyplot as plt

from plotter.findminmax import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--data', type=str, default="Q.dat", help='Name of timeseries to analyze.')
    parser.add_argument('--n_bins', type=int, default=20, help='Optional. number of bins for histogram.')
    parser.add_argument('--title', type=str, default="", help='Optional. Title for plot.')
    parser.add_argument('--saveas', type=str, default="None", help='Optional. Filename to save plot.')
    args = parser.parse_args()

    filename = args.data 
    coord_name = filename.split(".")[0]

    coord = np.loadtxt(filename)
    if filename == "Q.dat":
        coord += np.random.rand(coord.shape[0])

    n,bins = np.histogram(coord,bins=args.n_bins)
    F = -np.log(n)
    F -= min(F)
    mid_bin = 0.5*(bins[1:] + bins[:-1])

    x, F_fit, minxval, maxxval, omegamin, omegamax = from_timeseries(coord,bins=50,npoly=17)
    plt.plot(mid_bin,F,'o',color='b',label='data')
    plt.plot(x,F_fit(x),'b',lw=2,label='fit')

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

    for i in range(minxval.shape[0]):
        plt.axvline(x=min_state_bounds[i][0])
        plt.axvline(x=min_state_bounds[i][1])
    plt.legend()
    plt.savefig("testpmf.png")
    raise SystemExit



    plt.plot(x,F,'b',lw=2)
    plt.xlabel(coord_label,fontsize=16)
    plt.ylabel("F(%s) (k$_B$T)" % coord_label,fontsize=16)
    plt.title(args.title)
    if args.saveas != "None":
        plt.savefig(args.saveas)
    plt.show()

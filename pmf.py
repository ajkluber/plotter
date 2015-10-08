import argparse
import numpy as np
import matplotlib.pyplot as plt

import plotter.pmfutil as pmfutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--data', type=str, default="Q.dat", help='Name of timeseries to analyze.')
    parser.add_argument('--n_bins', type=int, default=60, help='Optional. number of bins for histogram.')
    parser.add_argument('--title', type=str, default="", help='Optional. Title for plot.')
    parser.add_argument('--saveas', type=str, help='Optional. Filename to save plot.')
    args = parser.parse_args()

    filename = args.data 
    coord_name = filename.split(".")[0]
    n_bins = args.n_bins
    title = args.title
    saveas = args.saveas

    coord = np.loadtxt(filename)
    if filename == "Q.dat":
        coord += np.random.normal(size=coord.shape[0])

    mid_bin, Fdata = pmfutil.pmf1D(coord,bins=n_bins)
    xinterp, F = pmfutil.interpolate_profile(mid_bin,Fdata)

    # Plot profile from data and interpolant 
    plt.plot(xinterp,F(xinterp),lw=2,color='b',label="data")
    plt.plot(mid_bin,Fdata,marker='o',color='b',label="fit")
    plt.legend()
    plt.xlabel("%s" % coord_name,fontsize=16)
    plt.ylabel("F(%s) (k$_B$T)" % coord_name,fontsize=16)
    plt.title(title)

    if saveas is not None:
        plt.savefig(saveas)
    plt.show()

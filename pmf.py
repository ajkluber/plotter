import argparse
import numpy as np
import matplotlib.pyplot as plt

#import plotter.pmfutil as pmfutil
import simulation.calc.pmfutil as pmfutil

def plotpmf1D(xvst,xlabel="",ylabel="Free energy (k$_B$T)",bins=50,saveas=None,display=True,label=""):
    """Plot 1D pmf"""
    if not display:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mid_bin,Fdata =  pmfutil.pmf1D(xvst,bins=bins) 
    plt.plot(mid_bin,Fdata,lw=2,label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if saveas is not None:
        plt.savefig(saveas) 
    if display:
        plt.show()

def plotpmf2D(xvst,yvst,xlabel="",ylabel="",bins=50,saveas=None,display=True):
    """Plot 1D pmf"""
    if not display:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xedges,yedges,Fdata = pmfutil.pmf2D(xvst,yvst,bins=bins) 
    plt.pcolormesh(xedges,yedges,Fdata)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    cbar = plt.colorbar()
    cbar.set_label("Free Energy (k$_B$T)")
    if saveas is not None:
        plt.savefig(saveas) 
    if display:
        plt.show()

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

    if filename.endswith("npy"):
        coord = np.load(filename)
    else:
        coord = np.loadtxt(filename)
    if filename == "Q.dat":
        coord += np.random.normal(size=coord.shape[0])

    mid_bin, Fdata = pmfutil.pmf1D(coord,bins=n_bins)
    #xinterp, F = pmfutil.interpolate_profile(mid_bin,Fdata)

    # Plot profile from data and interpolant 
    #plt.plot(xinterp,F(xinterp),lw=2,color='b',label="data")
    plt.plot(mid_bin,Fdata,marker='o',color='b',label="fit")
    plt.legend()
    plt.xlabel("%s" % coord_name,fontsize=16)
    plt.ylabel("F(%s) (k$_B$T)" % coord_name,fontsize=16)
    plt.title(title)

    if saveas is not None:
        plt.savefig(saveas)
    plt.show()

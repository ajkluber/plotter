import os
import argparse
import shutil
import logging
import numpy as np

import pmfutil

def save_Tf(abs_dFmin,uniq_Tlist,Tlist):
    """Save estimated folding temperature based on differen"""
    Tfidx = abs_dFmin.index(min(abs_dFmin))
    if abs_dFmin[Tfidx] > 2.0:
        print "difference in free energy is pretty large"
    Tf = uniq_Tlist[Tfidx]
    print "Estimating Tf as:", Tf
    with open("../Tf","w") as fout:
        fout.write(Tf)

    Tf_temps = ""
    for i in range(len(Tlist[Tfidx])):
        Tf_temps += "%s\n" % Tlist[Tfidx][i]
    with open("../Tf_temps","w") as fout:
        fout.write(Tf_temps)
    with open("Tf_temps","w") as fout:
        fout.write(Tf_temps)

    shutil.copy("%s_state_bounds.txt" % Tf,"state_bounds.txt")

def get_timeseries_by_temperature(temps_file,coord_file):
    """Get timeseries from directories, collected by unique temperatures"""
    coord_name = coord_file.split(".")[0]
    temps = [ x.rstrip("\n") for x in open(temps_file,"r").readlines() ]
    uniq_Tlist = []
    coordlist = []
    Tlist = []
    for i in range(len(temps)):
        T = temps[i].split("_")[-2]
        num = temps[i].split("_")[-1]
        if coord_file.endswith(".npy"):
            coordtemp = np.load("%s/%s" % (temps[i],coord_file))
        else:
            coordtemp = np.loadtxt("%s/%s" % (temps[i],coord_file))

        if coord_name == "Q":
            coordtemp += np.random.normal(size=len(coordtemp))

        if T not in uniq_Tlist:
            uniq_Tlist.append(T)
            coordlist.append(coordtemp)
            Tlist.append([temps[i]])
        else:
            idx = uniq_Tlist.index(T)
            coordlist[idx] = np.concatenate((coordlist[idx],coordtemp))
            Tlist[idx].append(temps[i])
    return uniq_Tlist, Tlist, coordlist

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--temps_file', type=str, required=True, help='File holding directory names.')
    parser.add_argument('--coord_file', type=str, default="Q.dat", help='Name of timeseries to analyze.')
    parser.add_argument('--n_bins', type=int, default=60, help='Optional. number of bins along coordinate.')
    parser.add_argument('--title', type=str, default="", help='Optional. Title for plot.')
    parser.add_argument('--nodisplay', action='store_true', help='Optional. Concatenate .')
    parser.add_argument('--saveas', type=str, default=None, help='Optional. Filename to save plot.')
    args = parser.parse_args()

    temps_file = args.temps_file 
    coord_file = args.coord_file 
    coord_name = coord_file.split(".")[0]
    coord_ext = coord_file.split(".")[-1]
    n_bins = args.n_bins 
    title = args.title
    saveas = args.saveas
    nodisplay = args.nodisplay

    if not os.path.exists("%s_profile" % coord_name):
        os.mkdir("%s_profile" % coord_name)

    # Get timeseries
    uniq_Tlist, Tlist, coordlist = get_timeseries_by_temperature(temps_file,coord_file)

    colors = ['b','r','g','k','cyan','magenta','salmon','darkgreen','gray']
    if nodisplay:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.chdir("%s_profile" % coord_name)

    plt.figure()
    abs_dFmin = []
    for i in range(len(coordlist)):
        T = uniq_Tlist[i]
        coordvst = coordlist[i]

        mid_bin, Fdata = pmfutil.pmf1D(coordvst,bins=n_bins)
        xinterp, F = pmfutil.interpolate_profile(mid_bin,Fdata)
        minidx, maxidx = pmfutil.extrema_from_profile(xinterp,F)
        min_bounds, max_bounds = pmfutil.state_bounds_from_profile(xinterp,F)
        min_labels, max_labels = pmfutil.assign_state_labels(min_bounds, max_bounds)
        pmfutil.save_state_bounds(T,min_bounds,max_bounds,min_labels,max_labels)

        dF_N_U = F(xinterp[minidx[-1]]) - F(xinterp[minidx[0]]) 
        abs_dFmin.append(abs(dF_N_U))

        # Determine highest barrier
        #barriers = [ F(xinterp[maxidx[i]]) - F(xinterp[minidx[0]]) for n in range(len(maxidx)) ]
        #max_barr = max(barriers) 
        #max_barr_idx = barriers.index(max_barr)
        #dF_TS_U = F(xinterp[minidx[max_barr_idx]]) - F(xinterp[minidx[0]]) 

        # Plot profile from data and interpolant 
        plt.plot(xinterp,F(xinterp),lw=2,color=colors[i])
        plt.plot(mid_bin,Fdata,'o',color=colors[i],label=T)

    save_Tf(abs_dFmin,uniq_Tlist,Tlist)

    plt.legend()
    plt.xlabel("%s" % coord_name,fontsize=16)
    plt.ylabel("F(%s) (k$_B$T)" % coord_name,fontsize=16)
    plt.title(title)

    if saveas is not None:
        plt.savefig(saveas)
    if not nodisplay:
        plt.show()

    os.chdir("..")


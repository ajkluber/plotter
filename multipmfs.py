import argparse
import numpy as np

#from bootFQ import get_F_with_error

def calculate_state_bounds(x,F,verbose=False):
    """Find the minima and maxima of inputted potential

        The bounds of a minimum end roughly where the free energy increases by
    ~kT to either side of the minimum.

    """

    xinterp = np.linspace(x.min(),x.max(),500)
    F_fit = np.poly1d(np.polyfit(x,F,9))
    dFdx = np.polyder(F_fit,m=1)
    A = dFdx(xinterp)
    minidx = np.where([(A[i] < 0) & (A[i + 1] > 0) for i in range(len(A) - 1) ])[0]
    maxidx = np.where([(A[i] > 0) & (A[i + 1] < 0) for i in range(len(A) - 1) ])[0]
    dF_minima = F_fit(xinterp[minidx[1]]) - F_fit(xinterp[minidx[0]])
    dF_barrier = F_fit(xinterp[maxidx[0]]) - F_fit(xinterp[minidx[0]])
    if len(minidx) != 2:
        print "WARNING! Found more/less than 2 minima!"
        print "    ",minidx
    state_bounds = []
    state_labels = ["U","N","TS"]
    
    for i in range(len(minidx)): 
        # Search to the left, then to the right around each minima
        left_min_bound = x.min()
        right_min_bound = x.max()
        for j in range(len(xinterp[:minidx[i]])):
            deltaF = F_fit(xinterp[minidx[i] - j]) - F_fit(xinterp[minidx[i]])
            if deltaF >= 0.3:
                left_min_bound = xinterp[minidx[i] - j]
                break
        for j in range(len(xinterp[minidx[i]:])):
            deltaF = F_fit(xinterp[minidx[i] + j]) - F_fit(xinterp[minidx[i]])
            if deltaF >= 0.3:
                right_min_bound = xinterp[minidx[i] + j]
                break
        state_bounds.append([left_min_bound,right_min_bound])
    left_max_bound = x.min()
    right_max_bound = x.max()
    for j in range(len(xinterp[:maxidx[0]])):
        deltaF = abs(F_fit(xinterp[maxidx[0] - j]) - F_fit(xinterp[maxidx[0]]))
        if deltaF >= 0.3:
            left_max_bound = xinterp[maxidx[0] - j]
            break
    for j in range(len(xinterp[maxidx[0]:])):
        deltaF = abs(F_fit(xinterp[maxidx[0] + j]) - F_fit(xinterp[maxidx[0]]))
        if deltaF >= 0.3:
            right_max_bound = xinterp[maxidx[0] + j]
            break
    state_bounds.append([left_max_bound,right_max_bound])

    if verbose:
        print "Found state bounds:"
        print "%8s   %8.2f %8.2f" % (state_labels[0],state_bounds[0][0],state_bounds[0][1])
        print "%8s   %8.2f %8.2f" % (state_labels[1],state_bounds[1][0],state_bounds[1][1])
        print "%8s   %8.2f %8.2f" % (state_labels[2],state_bounds[2][0],state_bounds[2][1])
    return xinterp, F_fit, dF_minima, dF_barrier, state_labels, state_bounds



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--coord_file', type=str, default="Q.dat", help='Name of timeseries to analyze.')
    parser.add_argument('--file', type=str, default="long_temps_last", help='File holding temps to plot.')
    parser.add_argument('--n_bins', type=int, default=25, help='Optional. number of bins along coordinate.')
    parser.add_argument('--n_histos', type=int, default=2, help='Optional. number of bootstrapping histograms.')
    parser.add_argument('--stride', type=int, default=1, help='Optional. number of frames to skip for subsampling.')
    parser.add_argument('--title', type=str, default="", help='Optional. Title for plot.')
    parser.add_argument('--all', action='store_true', help='Optional. Concatenate .')
    parser.add_argument('--verbose', action='store_true', help='Optional. Concatenate .')
    parser.add_argument('--nodisplay', action='store_true', help='Optional. Concatenate .')
    parser.add_argument('--saveas', type=str, default=None, help='Optional. Filename to save plot.')
    args = parser.parse_args()

    coord_file = args.coord_file 
    coord_name = coord_file.split(".")[0]
    coord_ext = coord_file.split(".")[-1]
    file = args.file
    n_bins = args.n_bins 
    n_histos = args.n_histos
    stride = args.stride
    title = args.title
    saveas = args.saveas
    verbose = args.verbose
    nodisplay = args.nodisplay

    if nodisplay:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = ['b','r','g','k','cyan','magenta','salmon','darkgreen']
        
    # Collect location of minima over temps
    # Determine state bounds. (1/3 kT)

    plt.figure()
    if args.all:
        # Sort temperatures
        temps = [ x.rstrip("\n") for x in open("long_temps","r").readlines() ]
        uniq_Tlist = []
        Qlist = []
        Tlist = []
        if verbose:
            print " Gathering unique Tlist"
        for i in range(len(temps)):
            T = temps[i].split("_")[0]
            num = temps[i].split("_")[1]
            if T not in uniq_Tlist:
                uniq_Tlist.append(T)
                Qtemp = np.loadtxt("%s/%s" % (temps[i],coord_file))
                if coord_name == "Q":
                    Qtemp += np.random.normal(size=len(Qtemp))
                Qlist.append(Qtemp)
                Tlist.append([temps[i]])
            else:
                idx = uniq_Tlist.index(T)
                Qtemp = np.loadtxt("%s/%s" % (temps[i],coord_file))
                if coord_name == "Q":
                    Qtemp += np.random.normal(size=len(Qtemp))
                Qlist[idx] = np.concatenate((Qlist[idx],Qtemp))
                Tlist[idx].append(temps[i])
        dFmin_all = []
        dFbarr_all = []
        state_bounds_all = []
        for i in range(len(Qlist)):
            # Subsample the coord_file. Calculate pmf.
            q_sub = Qlist[i][::stride]
            q_sub += np.random.normal(size=len(q_sub))
            n,bins = np.histogram(q_sub, bins=n_bins, density=False)
            bin_centers = 0.5*(bins[1:] + bins[:-1])
            F = -np.log(n)
            F -= min(F) 
            plt.plot(bin_centers,F,lw=2,color=colors[i],label=uniq_Tlist[i])

            xinterp, F_fit, dF_minima, dF_barrier, state_labels, state_bounds = calculate_state_bounds(bin_centers,F,verbose=verbose)
            state_bounds_all.append(state_bounds)
            dFmin_all.append(dF_minima)
            dFbarr_all.append(dF_barrier)
            plt.plot(xinterp,F_fit(xinterp),lw=2,color=colors[i])
            # Get minima, maxima, state_bounds.
        abs_dFmin_all = [ abs(x) for x in dFmin_all ]
        Tfidx = abs_dFmin_all.index(min(abs_dFmin_all))
        Tf = uniq_Tlist[Tfidx]
        if verbose:
            print "Estimating Tf as:", Tf
        Tf_temps = ""
        for i in range(len(Tlist[Tfidx])):
            Tf_temps += "%s\n" % Tlist[Tfidx][i]
        with open("Tf_temps","w") as fout:
            fout.write(Tf_temps)
        state_bounds = state_bounds_all[Tfidx]
        with open("state_bounds.txt","w") as fout:
            temp_string = "%s  %.2f %.2f\n" % (state_labels[0],state_bounds[0][0],state_bounds[0][1])
            temp_string += "%s  %.2f %.2f\n" % (state_labels[2],state_bounds[2][0],state_bounds[2][1])
            temp_string += "%s  %.2f %.2f\n" % (state_labels[1],state_bounds[1][0],state_bounds[1][1])
            fout.write(temp_string)
        if abs_dFmin_all[Tfidx] > 2.0:
            print "difference in free energy is pretty large"
    else:
        temps = [ x.rstrip("\n") for x in open(file,"r").readlines() ]
        for i in range(len(temps)):
            # Subsample the coord_file.
            q = np.loadtxt("%s/%s" % (temps[i],coord_file))
            if coord_name == "Q":
                q += np.random.normal(size=len(q))
            q_sub = q[::stride]

            # Calculate cumulative probability distribution of coordinate.
            n,bins = np.histogram(q_sub, bins=n_bins, density=False)
            bin_centers = 0.5*(bins[1:] + bins[:-1])
            F = -np.log(n)
            F -= min(F) 

            plt.plot(bin_centers,F,lw=2,color=colors[i],label=temps[i])

    plt.legend()
    plt.xlabel("%s" % coord_name,fontsize=20)
    plt.ylabel("F(%s) (k$_B$T)" % coord_name,fontsize=20)
    plt.title(title)
    if saveas is not None:
        plt.savefig(saveas)
    if not nodisplay:
        plt.show()

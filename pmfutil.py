import os
import numpy as np

def pmf1D(xvst,bins=50):
    """Histogram timeseries to get 1D pmf"""
    n,bins = np.histogram(xvst,bins=bins)
    mid_bin = 0.5*(bins[1:] + bins[:-1])
    Fdata = -np.log(n)
    Fdata -= min(Fdata)
    return mid_bin,Fdata

def pmf2D(xvst,yvst,bins=50):
    """Histogram timeseries to get 1D pmf"""
    nxy,xedges,yedges = np.histogram2d(xvst,yvst,bins=bins)
    Fdata = -np.log(nxy)
    Fdata -= Fdata.min()
    return xedges,yedges,Fdata

def interpolate_profile(mid_bin,Fdata,npoly=20,ninterp=500):
    """Interpolate 1D profile with polynomial"""
    xinterp = np.linspace(mid_bin.min(),mid_bin.max(),ninterp)
    F = np.poly1d(np.polyfit(mid_bin,Fdata,npoly))
    return xinterp, F

def interpolate_profiles(mid_bin,Fdatas,npoly=20,ninterp=500):
    """Interpolate 1D profile with polynomial"""
    xinterp = np.linspace(mid_bin.min(),mid_bin.max(),ninterp)
    Fs = []
    for i in range(Fdatas.shape[1]):
        Fs.append(np.poly1d(np.polyfit(mid_bin,Fdatas[:,i],npoly)))
    return xinterp, Fs

def cubic_interpolate_profiles(mid_bins,Fdatas,ninterp=500):
    """Interpolate 1D profile with polynomial"""
    from scipy.interpolate import interp1d
    Fs = []
    for i in range(len(Fdatas)):
        Fs.append(interp1d(mid_bins[i],Fdatas[i],kind="cubic"))
    return Fs


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
    # Determine state bounds for minima
    min_state_bounds = []
    for i in range(minidx.shape[0]):
        left_min_bound = xinterp.min()
        right_min_bound = xinterp.max()
        for j in range(xinterp[:minidx[i]].shape[0]):
            deltaF = F(xinterp[minidx[i] - j]) - F(xinterp[minidx[i]])
            if deltaF >= threshold:
                left_min_bound = xinterp[minidx[i] - j]
                break
        for j in range(xinterp[minidx[i]:].shape[0]):
            deltaF = F(xinterp[minidx[i] + j]) - F(xinterp[minidx[i]])
            if deltaF >= threshold:
                right_min_bound = xinterp[minidx[i] + j]
                break
        min_state_bounds.append([left_min_bound,right_min_bound])

    # Determine state bounds for maxima
    max_state_bounds = []
    for i in range(maxidx.shape[0]):
        left_min_bound = xinterp.min()
        right_min_bound = xinterp.max()
        for j in range(xinterp[:maxidx[i]].shape[0]):
            deltaF = abs(F(xinterp[maxidx[i] - j]) - F(xinterp[maxidx[i]]))
            if deltaF >= threshold:
                left_min_bound = xinterp[maxidx[i] - j]
                break
        for j in range(xinterp[maxidx[i]:].shape[0]):
            deltaF = abs(F(xinterp[maxidx[i] + j]) - F(xinterp[maxidx[i]]))
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
    counter = 1
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

def get_free_energy_profiles(sourceroot,parent_dirs,sub_dirs,coordfile,tempsfile,saveroot=None,nbins=40):
    """Get free energy profiles from coordinate"""
    cwd = os.getcwd()
    coordname = coordfile.split(".dat")[0]
    Fprofiles = [[] for i in range(len(parent_dirs)) ]
    Fmid_bins = [[] for i in range(len(parent_dirs)) ]
    for i in range(len(parent_dirs)):
        for j in range(len(sub_dirs)):
            os.chdir("%s/%s/%s" % (sourceroot,parent_dirs[i],sub_dirs[j]))
            if not os.path.exists("%s_profile/F.dat" % coordname):
            #if True:
                coordvst = np.concatenate([ np.loadtxt("%s/%s" % (x.rstrip("\n"),coordfile)) for x in open(tempsfile,"r").readlines() ])
                mid_bin, Fvsx = pmf1D(coordvst,bins=nbins)
                dFdx = np.array([ (Fvsx[x + 1] - Fvsx[x])/(mid_bin[x + 1] - mid_bin[x]) for x in range(len(Fvsx) - 1) ])
                minidx = np.where([(dFdx[x] < 0) & (dFdx[x + 1] > 0) for x in range(len(dFdx) - 1) ])[0]
                maxidx = np.where([(dFdx[x] > 0) & (dFdx[x + 1] < 0) for x in range(len(dFdx) - 1) ])[0]
                
                if saveroot is not None:
                    os.chdir("%s/%s/%s" % (saveroot,parent_dirs[i],sub_dirs[j]))
                if not os.path.exists("%s_profile" % coordname):
                    os.mkdir("%s_profile" % coordname)
                os.chdir("%s_profile" % coordname)
                np.savetxt("F.dat",Fvsx)
                np.savetxt("mid_bin.dat",mid_bin)
                np.savetxt("maxima.dat",mid_bin[maxidx])
                np.savetxt("minima.dat",mid_bin[minidx])
            else:
                Fvsx = np.loadtxt("%s_profile/F.dat" % coordname)
                mid_bin = np.loadtxt("%s_profile/mid_bin.dat" % coordname)
                dFdx = np.array([ (Fvsx[x + 1] - Fvsx[x])/(mid_bin[x + 1] - mid_bin[x]) for x in range(len(Fvsx) - 1) ])
                minidx = np.where([(dFdx[x] < 0) & (dFdx[x + 1] > 0) for x in range(len(dFdx) - 1) ])[0]
                maxidx = np.where([(dFdx[x] > 0) & (dFdx[x + 1] < 0) for x in range(len(dFdx) - 1) ])[0]

                if saveroot is not None:
                    os.chdir("%s/%s/%s" % (saveroot,parent_dirs[i],sub_dirs[j]))
                    if not os.path.exists("%s_profile" % coordname):
                        os.mkdir("%s_profile" % coordname)
                    os.chdir("%s_profile" % coordname)
                    np.savetxt("F.dat",Fvsx)
                    np.savetxt("mid_bin.dat",mid_bin)
                    np.savetxt("maxima.dat",mid_bin[maxidx])
                    np.savetxt("minima.dat",mid_bin[minidx])

            Fprofiles[i].append(Fvsx)
            Fmid_bins[i].append(mid_bin)
    os.chdir(cwd)
    return Fmid_bins, Fprofiles

def gridplot_Fvsx(myroot,parent_dirs,sub_dirs,Fmid_bins,Fprofiles,Fmax,coordname,coordsymb,name,display=True):
    """Plot free energy profiles in grid"""
    if not display:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dFdagg = [[] for i in range(len(parent_dirs)) ]
    for i in range(len(parent_dirs)):
        counter = 0
        fig1,axes = plt.subplots(3,4,sharex=True,sharey=True)
        fig2,ax2 = plt.subplots(1,1)
        for j in range(len(sub_dirs)):
            Fvsx = Fprofiles[i][j]
            mid_bin = Fmid_bins[i][j]
            ax2.plot(mid_bin,Fvsx,label="replica %d" % (j+1))

            # Determine minima locations
            dFdx = np.array([ (Fvsx[x + 1] - Fvsx[x])/(mid_bin[x + 1] - mid_bin[x]) for x in range(len(Fvsx) - 1) ])
            minidx = np.where([(dFdx[x] < 0) & (dFdx[x + 1] > 0) for x in range(len(dFdx) - 1) ])[0] + 1
            maxidx = np.where([(dFdx[x] > 0) & (dFdx[x + 1] < 0) for x in range(len(dFdx) - 1) ])[0] + 1
            dFdagg[i].append(Fvsx[maxidx[0]] - Fvsx[minidx[0]])

            i_idx = counter / 4
            j_idx = counter % 4
            ax = axes[i_idx,j_idx]

            ax.plot(mid_bin,Fvsx,color="#5DA5DA")
            ax.axvline(mid_bin[minidx[0]],ymin=0,ymax=0.2)
            ax.axvline(mid_bin[minidx[-1]],ymin=0,ymax=0.2)
            ax.axvline(mid_bin[maxidx[0]],ymin=0,ymax=0.2)
            ax.set_ylim(0,Fmax)
            
            counter += 1
        axes[1,0].set_ylabel("$F(%s)$ (k$_B$T)" % coordsymb.replace("$",""))
        axes[2,1].set_xlabel(coordsymb)
        fig1.suptitle("Free energy curves %s" % name,fontsize=18)
        fig1.subplots_adjust(hspace=0,wspace=0)

        ax2.legend()
        ax2.set_title("Replica free energy profiles %s" % name)
        ax2.set_xlabel(coordsymb)
        ax2.set_ylabel("$F(%s)$ (k$_B$T)" % coordsymb.replace("$",""))
        if not os.path.exists("plots"):
            os.mkdir("plots")
        fig2.savefig("plots/all_F_vs_%s.png" % coordname)
        fig1.savefig("plots/all_F_vs_%s_grid.png" % coordname)

        np.savetxt("%s/%s/%s_dFdagg.dat" % (myroot,parent_dirs[i],coordname),np.array(dFdagg[i]))
    return dFdagg

def plot_dFdagg_vs_B(dFdagg,Fmax,variance,name,coordname,display=True):
    """Plot free energy barrier heights versus b """
    if not display:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.figure()
    for i in range(len(dFdagg)):
        plt.errorbar(np.sqrt(float(variance[i])),np.mean(dFdagg[i]),yerr=np.std(dFdagg[i]),color="#5DA5DA")
        plt.plot(np.sqrt(float(variance[i])),np.mean(dFdagg[i]),marker='o',color="#5DA5DA")
    plt.ylabel("Folding free energy barrier (k$_B$T)")
    plt.xlabel("Frustration $b$")
    plt.title("Barrier Heights %s" % name)
    plt.ylim(0,Fmax)
    plt.savefig("plots/%s_dFdagg_vs_b.png" % coordname)

def plot_Fvsx_variance_gridplot(myroot,parent_dirs,sub_dirs,variance,Fprofiles,Fmid_bins,coloridx,name,coordname,display=True):
    """Gridplot free energy profiles with average profile in bold """
    if not display:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plotter.cube_cmap import cubecmap
    from scipy.interpolate import interp1d
    fig1,axes = plt.subplots(3,4,sharex=True,sharey=True)
    counter = 0
    Finterp_all = []
    for i in range(len(parent_dirs)):
        Finterp_rep = []
        xall = []
        i_idx = counter / 4
        j_idx = counter % 4
        ax = axes[i_idx,j_idx]
        for j in range(len(sub_dirs)):
            Fvsx = Fprofiles[i][j]
            mid_bin = Fmid_bins[i][j]

            ax.plot(mid_bin,Fvsx,color=cubecmap(coloridx[i]),alpha=0.2)

            Finterp = interp1d(mid_bin,Fvsx,kind="cubic")
            Finterp_rep.append(Finterp)
            xall.append(mid_bin)

        Finterp_all.append(Finterp_rep)
        ax.text(mid_bin.max()*0.3,3.2,"$b^2 = %s$" % variance[i])
        qmax = min([ max(xall[x]) for x in range(len(xall)) ])
        qmin = max([ min(xall[x]) for x in range(len(xall)) ])
        q = np.linspace(qmin,qmax,100)
        Favg = np.mean(np.array([ Finterp_rep[x](q) for x in range(len(Finterp_rep)) ]),axis=0)
        #ax.plot(q,Favg,color="#5DA5DA")
        ax.plot(q,Favg,color=cubecmap(coloridx[i]),lw=3)
        ax.set_ylim(0,4)

        counter += 1

    axes[1,0].set_ylabel("Free energy F(Q) (k$_B$T)")
    axes[2,1].set_xlabel("Folding Q")
    fig1.suptitle("Free energy curves %s" % name,fontsize=18)
    fig1.subplots_adjust(hspace=0,wspace=0)
    if not os.path.exists("plots"):
        os.mkdir("plots")
    fig1.savefig("plots/all_b2_profiles_%s.png" % coordname)

def calculate_dFdagg(myroot,parent_dirs,sub_dirs,Fmid_bins,Fprofiles,coordname):
    """Calculate free energy barrier height"""
    from scipy.interpolate import interp1d
    dFdagg = [[] for i in range(len(parent_dirs)) ]
    Finterp = [[] for i in range(len(parent_dirs)) ]
    for i in range(len(parent_dirs)):
        for j in range(len(sub_dirs)):
            Fvsx = Fprofiles[i][j]
            mid_bin = Fmid_bins[i][j]

            F = interp1d(mid_bin,Fvsx,kind="cubic")
            Finterp.append(F)

            # Determine minima locations
            dFdx = np.array([ (Fvsx[x + 1] - Fvsx[x])/(mid_bin[x + 1] - mid_bin[x]) for x in range(len(Fvsx) - 1) ])
            minidx = np.where([(dFdx[x] < 0) & (dFdx[x + 1] > 0) for x in range(len(dFdx) - 1) ])[0] + 1
            maxidx = np.where([(dFdx[x] > 0) & (dFdx[x + 1] < 0) for x in range(len(dFdx) - 1) ])[0] + 1
            dFdagg[i].append(Fvsx[maxidx[0]] - Fvsx[minidx[0]])

        np.savetxt("%s/%s/%s_dFdagg.dat" % (myroot,parent_dirs[i],coordname),np.array(dFdagg[i]))
    return dFdagg, Finterp

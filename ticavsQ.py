import os
import argparse
import numpy as np 
import matplotlib.pyplot as plt

from plotter.cube_cmap import cubecmap
from plotter.shiftcmap import shiftedColorMap

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def misc():
    # Q vs TICA 1
    z,x,y = np.histogram2d(Qw, tica, bins=(90,90))
    F = -np.log(z)
    F -= F.min()
    extent = [x[0], x[-1], y[0], y[-1]]
    ax1 = plt.contourf(F.T, 100, extent=extent, cmap=cubecmap)
    xlabel("Q",fontsize=16)
    ylabel("TICA 1",fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Free energy')


    # TICA1 vs TICA2 (if possible)
    plt.figure()
    z,x,y = np.histogram2d(Q, tica, bins=(90,90))
    F = -np.log(z)
    F -= F.min()
    extent = [x[0], x[-1], y[0], y[-1]]
    ax1 = plt.contourf(F.T, 100, extent=extent)
    xlabel("Q",fontsize=16)
    ylabel("TICA 1",fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Free energy')

    
    plt.figure()
    z,x,y = np.histogram2d(tica, newtica, bins=(90,90))
    F = -np.log(z)
    F -= F.min()
    extent = [x[0], x[-1], y[0], y[-1]]
    ax1 = plt.contourf(F.T, 100, extent=extent)
    xlabel("Q",fontsize=16)
    ylabel("TICA 1",fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Free energy')

    
    figure()
    z,x,y = histogram2d(tica1, tica2, bins=(90,90))
    F = -log(z)
    F -= F.min()
    extent = [x[0], x[-1], y[0], y[-1]]
    ax1 = contourf(F.T, 100, extent=extent, cmap=cubecmap)
    xlabel("TICA 1",fontsize=16)
    ylabel("TICA 2",fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Free energy')

    figure()
    z,x,y = histogram2d(Q, tica1, bins=(90,90))
    F = -log(z)
    F -= F.min()
    extent = [x[0], x[-1], y[0], y[-1]]
    ax1 = contourf(F.T, 100, extent=extent, cmap=cubecmap)
    xlabel("Q",fontsize=16)
    ylabel("TICA 2",fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Free energy')

def plot_tica_eigenvalues(s,lagtau,saveas=True):
    first_neg = np.where(s < 0)[0][0]
    imp_timescales = -lagtau/np.log(s[:first_neg])

    # Plot eigenvalues
    plt.figure()
    plt.plot(s,'g.')
    plt.axvline(x=first_neg,color='k')
    plt.title("First negative %d" % first_neg,fontsize=16)
    plt.xlabel("Eigenvalue index",fontsize=16)
    plt.ylabel("Eigenvalue",fontsize=16)
    if saveas:
        plt.savefig("eigenvalues.png",bbox_inches="tight")
        plt.savefig("eigenvalues.pdf",bbox_inches="tight")
        plt.savefig("eigenvalues.eps",bbox_inches="tight")

    # Implied timescales
    plt.figure()
    plt.plot(imp_timescales,'g.')
    plt.title("Implied timescales")
    plt.xlabel("Eigenvalue index $i$",fontsize=16)
    plt.ylabel("Implied timescale  $t_i$ (frames)",fontsize=16)
    if saveas:
        plt.savefig("implied_timescales.png",bbox_inches="tight")
        plt.savefig("implied_timescales.pdf",bbox_inches="tight")
        plt.savefig("implied_timescales.eps",bbox_inches="tight")
    
    # Ratio of subsequent timescales
    plt.figure()
    plt.plot(imp_timescales[:-1]/imp_timescales[1:],'g.')
    plt.title("Ratio of subsequent implied timescales",fontsize=16)
    plt.xlabel("Eigenvalue index $i$",fontsize=16)
    plt.ylabel("Implied timescale $t_i/t_{i+1}$",fontsize=16)
    if saveas:
        plt.savefig("ratio_implied_timescales.png",bbox_inches="tight")
        plt.savefig("ratio_implied_timescales.pdf",bbox_inches="tight")
        plt.savefig("ratio_implied_timescales.eps",bbox_inches="tight")
    plt.show()


def plot_2Dpmf(X,Y,xlabel,ylabel,title,saveas,bins=(100,100)):
    cubecmap.set_bad(color="lightgray",alpha=1.)
    plt.figure()
    z,x,y = np.histogram2d(X,Y,bins=bins)
    F = -np.log(z)
    F -= F.min()
    extent = [x[0], x[-1], y[0], y[-1]]
    F = np.ma.masked_invalid(F)
    ax1 = plt.contourf(F.T,100,extent=extent,cmap=cubecmap)
    plt.xlabel(xlabel,fontsize=16)
    plt.ylabel(ylabel,fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Free energy (k$_B$T)',fontsize=15)
    plt.title(title,fontsize=16)
    if saveas is not None:
        plt.savefig(saveas+".png",bbox_inches="tight")
        plt.savefig(saveas+".pdf",bbox_inches="tight")
        plt.savefig(saveas+".eps",bbox_inches="tight")

def plot_1Dpmf(X,xlabel,title,saveas,nbins=100):
    n,bins = np.histogram(X,bins=nbins)
    pmf = -np.log(n)
    pmf -= pmf.min()
    plt.figure()
    plt.plot(0.5*(bins[1:] + bins[:-1]),pmf,lw=2)
    plt.xlabel(xlabel,fontsize=16)
    plt.ylabel("Free energy (k$_B$T)",fontsize=16)
    plt.title(title,fontsize=16)
    if saveas is not None:
        plt.savefig(saveas+".png",bbox_inches="tight")
        plt.savefig(saveas+".pdf",bbox_inches="tight")
        plt.savefig(saveas+".eps",bbox_inches="tight")

def plot_two_1Dpmfs(X1,X2,xlabel,X1label,X2label,title,saveas,nbins=100):
    X1 -= X1.min()
    X1 /= X1.max()
    plt.figure()
    X2 /= X2.max()
    n2,bins2 = np.histogram(X2,bins=nbins)
    pmf2 = -np.log(n2)
    pmf2 -= pmf2.min()
    plt.plot(0.5*(bins2[1:] + bins2[:-1]),pmf2,'r',lw=2,label=X2label)

    n1,bins1 = np.histogram(X1,bins=nbins)
    pmf1 = -np.log(n1)
    pmf1 -= pmf1.min()
    plt.plot(0.5*(bins1[1:] + bins1[:-1]),pmf1,'b',lw=2,label=X1label)
    plt.xlabel(xlabel,fontsize=16)
    plt.ylabel("Free energy (k$_B$T)",fontsize=16)
    plt.title(title,fontsize=16)
    plt.legend(prop={'size':16})
    if saveas is not None:
        plt.savefig(saveas+".png",bbox_inches="tight")
        plt.savefig(saveas+".pdf",bbox_inches="tight")
        plt.savefig(saveas+".eps",bbox_inches="tight")

def plot_weights(C,coordname,title,saveas,shift):
    C = np.ma.masked_invalid(C)
    cubecmap.set_bad(color="lightgray",alpha=1.)
    rwb = plt.get_cmap("bwr_r")
    shifted_cmap = shiftedColorMap(rwb, midpoint=shift, name='shifted')
    plt.figure()
    #plt.pcolormesh(C,cmap=cubecmap)
    plt.pcolormesh(C,cmap=shifted_cmap)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Contact weight on %s' % coordname,fontsize=15)
    plt.xlabel("Residue i",fontsize=16)
    plt.ylabel("Residue j",fontsize=16)
    plt.xlim(0,C.shape[0])
    plt.ylim(0,C.shape[0])
    plt.title(title,fontsize=16)
    if saveas is not None:
        plt.savefig(saveas+".png",bbox_inches="tight")
        plt.savefig(saveas+".pdf",bbox_inches="tight")
        plt.savefig(saveas+".eps",bbox_inches="tight")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--temps', type=str, required=True, help='File with temps to use.')
    parser.add_argument('--feature', type=str, required=True, help='all_contacts or native_contacts.')
    parser.add_argument('--lag', type=int, required=True, help='Lag for TICA.')
    parser.add_argument('--stride', type=int, required=True, help='Name of directory.')
    parser.add_argument('--n_bins', type=int, default=100, help='Name of directory.')
    parser.add_argument('--no_display', action='store_true', help="Don't show plots.")
    args = parser.parse_args()

    temps = args.temps
    feature = args.feature
    lag = args.lag
    stride = args.stride
    no_display = args.no_display

    dirs = [ x.rstrip("\n") for x in  open(temps,"r").readlines() ]
    n_residues = len(open("%s/Native.pdb" % dirs[0],"r").readlines()) - 1

    available_features = ["native_contacts","all_contacts"]
    if feature not in available_features:
        raise IOError("--feature should be in: %s" % available_features.__str__())

    if feature == "all_contacts":
        prefix = "all"
    else:
        prefix = "nat"

    # Get tica coordinates  
    psi1 = np.concatenate([ np.loadtxt("%s/tica1_%s_%d_%d.dat" % (x,prefix,lag,stride)) for x in dirs ])
    if os.path.exists("%s/tica2_%s_%d_%d.dat" % (dirs[0],prefix,lag,stride)):
        psi2 = np.concatenate([ np.loadtxt("%s/tica2_%s_%d_%d.dat" % (x,prefix,lag,stride)) for x in dirs ])
    else:
        psi2 = None
    coord = np.concatenate([ np.loadtxt("%s/Q.dat" % x) + np.random.rand(file_len("%s/Q.dat" % x)) for x in dirs ])
    corr = np.dot(coord,psi1)/(np.linalg.norm(coord)*np.linalg.norm(psi1))

    # Scale tica1 coord so it is (0,1) and correlated with Q.
    A = np.sign(corr)*psi1.max()
    psi1 /= A

    os.chdir("tica_%s_%d_%d" % (prefix,lag,stride))
    with open("Qcorrsign","w") as fout:
        fout.write("%d" % np.sign(corr))

    # Get tica weights
    psi1_w = np.loadtxt("tica1_%s_%d_%d_weights.dat" % (prefix,lag,stride))
    psi1_w[:,2] /= A
    psi1 -= psi1.min()
    psi1 /= psi1.max()
    Cpsi1 = np.zeros((n_residues,n_residues))*np.nan
    for i in range(psi1_w.shape[0]):
        Cpsi1[int(psi1_w[i,1]),int(psi1_w[i,0])] = psi1_w[i,2]

    if psi2 is not None:
        psi2_w = np.loadtxt("tica2_%s_%d_%d_weights.dat" % (prefix,lag,stride))
        Cpsi2 = np.zeros((n_residues,n_residues))*np.nan
        for i in range(psi2_w.shape[0]):
            Cpsi2[int(psi2_w[i,1]),int(psi2_w[i,0])] = psi2_w[i,2]

    tica_eigenvalues = np.loadtxt("eigenvalues.dat")

    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")

    # Plot tica eigenvalues and implied timescales.
    plot_tica_eigenvalues(tica_eigenvalues,lag)

    # Plot 1D 
    plot_1Dpmf(coord,"Native contacts $Q$","Native contacts","Qpmf")
    plot_1Dpmf(psi1,"Foldedness  $\widetilde{\psi}_1$","TICA 1 $\widetilde{\psi}_1$","psi1pmf")
    plot_two_1Dpmfs(psi1,coord,"Foldedness","$\widetilde{\psi}_1$","$Q$","$F(Q)$ and $F(\widetilde{\psi}_1)$","Qandpsi1pmf")

    shift1 = 1 - (psi1_w[:,2].max()/(psi1_w[:,2].max() + abs(psi1_w[:,2].min())))
    plot_weights(Cpsi1,"$\widetilde{\psi}_1$","Contact weights for $\widetilde{\psi}_1$","psi1_w",shift1)

    # Plot tica1 versus Q
    plot_2Dpmf(coord,psi1,"$Q$","TICA 1 $\widetilde{\psi}_1$","TICA 1 vs Q","Qvspsi1")

    # Plot tica2 vs tica1, tica2 vs Q
    if psi2 is not None:
        shift2 = 1 - (psi2_w[:,2].max()/(psi2_w[:,2].max() + abs(psi2_w[:,2].min())))
        plot_weights(Cpsi2,"$\widetilde{\psi}_2$","Contact weights for $\widetilde{\psi}_2$","psi2_w",shift2)
        plot_1Dpmf(psi2,"Second relaxtion $\widetilde{\psi}_2$","TICA 2 $\widetilde{\psi}_2$","psi2pmf")
        plot_2Dpmf(coord,psi2,"$Q$","TICA 2 $\widetilde{\psi}_2$","TICA 2 vs Q","Qvspsi2")
        plot_2Dpmf(psi1,psi2,"TICA 1 $\widetilde{\psi}_1$","TICA 2 $\widetilde{\psi}_2$","TICA 1 vs TICA 2","psi1vspsi2")

    if not no_display:
        plt.show()

    os.chdir("../..")

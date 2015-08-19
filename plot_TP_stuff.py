import os
import numpy as np
import matplotlib.pyplot as plt

from misc.analysis.TPT import partition_TP
from misc.acf import calculate_acf

def plot_TP_quantitites(dwellA,dwellB,TP_A_B,TP_B_A):
    dwellA_times = [ len(A) for A in dwellA ]
    dwellB_times = [ len(A) for A in dwellB ]
    TP_A_B_times = [ len(A) for A in TP_A_B ]
    TP_B_A_times = [ len(A) for A in TP_B_A ]
    avgTP_A_B_time = np.mean(TP_A_B_times)
    avgTP_B_A_time = np.mean(TP_B_A_times)
    avgTP_time = np.mean(TP_A_B_times + TP_B_A_times)

    TP_A_B_time_dist, TP_A_Bbins = np.histogram(TP_A_B_times,bins=30)
    TP_A_B_bin_centers = 0.5*(TP_A_Bbins[1:] + TP_A_Bbins[:-1])

    TP_B_A_time_dist, TP_B_Abins = np.histogram(TP_B_A_times,bins=30)
    TP_B_A_bin_centers = 0.5*(TP_B_Abins[1:] + TP_B_Abins[:-1])

    TP_time_dist, TPbins = np.histogram(TP_A_B_times + TP_B_A_times,bins=30)
    TPbin_centers = 0.5*(TPbins[1:] + TPbins[:-1])

    N = float(len(x))
    N_TP = float(sum(TP_A_B_times) + sum(TP_B_A_times))

    Nr, bins = np.histogram(x,bins=bins)
    Nr_TP, bins = np.histogram(np.concatenate(TP_A_B + TP_B_A),bins=bins)
    bin_centers = 0.5*(bins[1:] + bins[:-1])

    Pr = Nr/N
    Pr_TP = Nr_TP/N_TP
    P_TP = N_TP/N
    P_TP_r = Pr_TP*P_TP/Pr

    pmf = -np.log(Pr)
    pmf -= pmf.min()

    plt.figure()
    plt.fill_between(bin_centers,pmf,facecolor='b',alpha=0.5)
    plt.plot(bin_centers,pmf,lw=2)
    plt.xlabel("Reaction Coordinate x",fontsize=16)
    plt.ylabel("Potential mean force F(x) (k$_B$T)",fontsize=16)
    plt.title("Toy model double well potential",fontsize=14)
    plt.xlim(min(bin_centers),max(bin_centers))
    plt.savefig("pmf_x.png",bbox_inches="tight")
    plt.savefig("pmf_x.pdf",bbox_inches="tight")
    plt.savefig("pmf_x.eps",bbox_inches="tight")


    plt.figure()
    plt.fill_between(bin_centers,P_TP_r,facecolor='b',alpha=0.5)
    plt.plot(bin_centers,P_TP_r,lw=2)
    plt.xlabel("Reaction Coordinate x",fontsize=16)
    plt.ylabel("P(TP|x)",fontsize=16)
    plt.title("Good reaction coordinates have one peak",fontsize=14)
    plt.savefig("P_TP_x.png",bbox_inches="tight")
    plt.savefig("P_TP_x.pdf",bbox_inches="tight")
    plt.savefig("P_TP_x.eps",bbox_inches="tight")

    plt.figure()
    plt.fill_between(bin_centers,Pr,facecolor='b',alpha=0.5)
    plt.plot(bin_centers,Pr,lw=2,color='b')
    plt.xlabel("Reaction Coordinate x",fontsize=16)
    plt.ylabel("P(x)",fontsize=16)
    plt.savefig("Px.png",bbox_inches="tight")
    plt.savefig("Px.pdf",bbox_inches="tight")
    plt.savefig("Px.eps",bbox_inches="tight")

    plt.figure()
    plt.fill_between(TPbin_centers,TP_time_dist,facecolor='b',alpha=0.5)
    plt.plot(TPbin_centers,TP_time_dist,'b',lw=2)
    plt.axvline(x=avgTP_time,color='k')
    plt.xlabel("Transition path time $t_{TP}$ (frames)",fontsize=16)
    plt.ylabel("P($t_{TP}$)",fontsize=16)
    plt.title("$\\langle t_{TP} \\rangle = %.2f$ frames " % avgTP_time,fontsize=16)
    plt.savefig("TP_time_dist.png",bbox_inches="tight")
    plt.savefig("TP_time_dist.pdf",bbox_inches="tight")
    plt.savefig("TP_time_dist.eps",bbox_inches="tight")

    plt.figure()
    plt.fill_between(TP_A_B_bin_centers,TP_A_B_time_dist,facecolor='b',alpha=0.5)
    plt.plot(TP_A_B_bin_centers,TP_A_B_time_dist,'b',lw=2)
    plt.axvline(x=avgTP_time,color='k')
    plt.xlabel("Transition path time $t_{TP}$ (frames)",fontsize=16)
    plt.ylabel("P($t_{TP}$)",fontsize=16)
    plt.title("Forward transition paths $\\langle t_{TP} \\rangle = %.2f $ " % avgTP_A_B_time)
    plt.savefig("TP_A_B_time_dist.png",bbox_inches="tight")
    plt.savefig("TP_A_B_time_dist.pdf",bbox_inches="tight")
    plt.savefig("TP_A_B_time_dist.eps",bbox_inches="tight")

    plt.figure()
    plt.fill_between(TP_B_A_bin_centers,TP_B_A_time_dist,facecolor='b',alpha=0.5)
    plt.plot(TP_B_A_bin_centers,TP_B_A_time_dist,'b',lw=2)
    plt.axvline(x=avgTP_time,color='k')
    plt.xlabel("Transition path time $t_{TP}$ (frames)",fontsize=16)
    plt.ylabel("P($t_{TP}$)",fontsize=16)
    plt.title("Backward transition paths $\\langle t_{TP} \\rangle = %.2f $ " % avgTP_B_A_time)
    plt.savefig("TP_B_A_time_dist.png",bbox_inches="tight")
    plt.savefig("TP_B_A_time_dist.pdf",bbox_inches="tight")
    plt.savefig("TP_B_A_time_dist.eps",bbox_inches="tight")

if __name__ == "__main__":

    nbins = 40

    os.chdir("double_well_1")

    print "loading data"
    x = np.load("x_1.npy")
    #v = np.load("v_0.npy")

    #acf = calculate_acf(x)
    # 1E7 frames
    #bins = np.linspace(x.min(),x.max(),nbins)
    bins = np.linspace(-8,8,nbins)

    #stateA = -4.6012024
    #stateB = 5.40280561
    stateA = -5
    stateB = 5

    print "paritioning trajectory into transition paths"
    dwellA, dwellB, TP_A_B, TP_B_A = partition_TP(x,stateA,stateB)

    plot_TP_quantitites(dwellA,dwellB,TP_A_B,TP_B_A)

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

import plotter.pmfutil as pmfutil
from plotter.cube_cmap import cubecmap
import simulation.calc.binned_contacts 

""" Mechanism as a function of heterogeneity

Heterogeneity should perturb the homogeneous mechanism.

At what value of heterogeneity do given topologies .

Folding mechanism should broaden from the homogeneous mechanism. perhaps

"""

colors2 = [('lightgray')] + [(cm.Blues(i)) for i in xrange(1,256)]
global gray_background
gray_background = matplotlib.colors.LinearSegmentedColormap.from_list('gray_background', colors2, N=256)

class ContactArgs(object):
    """ Data object to parameterize contact function """
    def __init__(self,trajsfile,bins=40,function="tanh",tanh_scale=0.05,
            coordfile="Qtanh_0_05.dat",topology="Native.pdb",chunksize=1000,periodic=False):
        self.trajs = trajsfile
        self.function = function
        self.coordfile = coordfile
        self.bins = bins
        self.tanh_scale = tanh_scale
        self.chunksize = chunksize
        self.topology = topology
        self.periodic = periodic

def plot_replica_TS_grid(parent_dirs,sub_dirs,pairs,n_residues,coordfile,tempsfile,tempsfile_backup,save_formats=["png"]):
    """ Plot grid of TS contact probabilities for replicas


    """
    coordname = coordfile.split(".")[0]
    n_pairs = pairs.shape[0]
    top_cwd = os.getcwd()
    all_rep_avg = []
    for i in range(len(parent_dirs)):
        # Plot all replicas for a given parent directory
        parent_dir = parent_dirs[i]
        os.chdir(parent_dir)
        parent_cwd = os.getcwd()
        print parent_dir

        fig1,axes = plt.subplots(3,4,sharex=True,sharey=True)

        counter = 0
        for rep_idx in range(len(sub_dirs)):
            # Plot one replica
            sub_dir = sub_dirs[rep_idx]
            os.chdir(sub_dir)

            i_idx = counter / 4
            j_idx = counter % 4
            ax = axes[i_idx,j_idx]

            ax.text(n_residues*0.4,n_residues*0.1,"replica %d" % (rep_idx + 1),fontsize=10)
            TS = plot_single_replica_TS(coordfile,tempsfile,tempsfile_backup)

            if TS is not None: 
                C = np.zeros((n_residues,n_residues))
                for j in range(pairs.shape[0]):
                    C[pairs[j,1],pairs[j,0]] = TS[j]
                # Plot TS contact map
                ax.pcolormesh(C,vmin=0,vmax=1.,cmap=gray_background)
                ax.plot([0,n_residues],[0,n_residues],'k',lw=1)
                ax.grid(True)
                ax.set_xlim(0,n_residues)
                ax.set_ylim(0,n_residues)

            counter += 1
            os.chdir(parent_cwd)

        if not os.path.exists("plots"):
            os.mkdir("plots")
        os.chdir("plots")
        # Save grid figure of TS's for all replicas
        fig1.suptitle("TS contact probability %s" % variance[i],fontsize=18)
        fig1.subplots_adjust(hspace=0,wspace=0)
        for format in save_formats:
            fig1.savefig("all_rep_TS.%s" % format,bbox_inches="tight")
        plt.close(fig1)

        os.chdir("..")
        os.chdir(top_cwd)

def plot_single_replica_TS(coordfile,tempsfile,tempsfile_backup):
    """ Plot single replica TS contact probabilities
    
    
    """
    coordname = coordfile.split(".")[0]

    # Check that reaction coordinate exists.
    if os.path.exists(tempsfile):
        with open(tempsfile,"r") as fin:
            Tdirs = [ x.rstrip("\n") for x in fin.readlines() ]
        coordfiles_exist = all([ os.path.exists("%s/%s" % (x,coordfile)) for x in Tdirs ])
    elif os.path.exists(tempsfile_backup):
        with open(tempsfile_backup,"r") as fin:
            Tdirs = [ x.rstrip("\n") for x in fin.readlines() ]
        coordfiles_exist = all([ os.path.exists("%s/%s" % (x,coordfile)) for x in Tdirs ])
    else:
        coordfiles_exist = False

    if not coordfiles_exist: 
        print "replica doesn't have %s" % coordfile
        TS = None
    else:
        # Load or calculate TS contact probabilities
        if os.path.exists("binned_contacts_vs_%s/cont_prob_TS.dat" % coordname):
            print "loading TS"
            TS = np.loadtxt("binned_contacts_vs_%s/cont_prob_TS.dat" % coordname)
        else:
            contact_args = ContactArgs("Tftrajs")
            TS = calculate_TS_probabilities(Tdirs,coordfile,contact_args)
    return TS

def calculate_TS_probabilities(Tdirs,coordfile,contact_args):
    """ Calculate the TS contact probabilities along some reaction coordinate

    Parameters
    ----------

    Tdirs : list
        List of directory names to collect reaction coordinate file from

    coordfile : str
        Name of file containing reaction coordinate.

    contact_args : object
        Object containing the 

    """

    T = Tdirs[0].split("_")[0]
    coordname = coordfile.split(".")[0]
    if not os.path.exists("%s_profile/state_bounds.txt" % coordname):
        # Determine state bounds from 1D profile
        print "calculating state bounds"
        coordvst = np.concatenate([np.loadtxt("%s/%s" % (x,coordfile)) for x in Tdirs ])
        os.chdir("%s_profile" % coordname)
        mid_bin, Fdata = pmfutil.pmf1D(coordvst,bins=40)
        xinterp, F = pmfutil.interpolate_profile(mid_bin,Fdata)
        minidx, maxidx = pmfutil.extrema_from_profile(xinterp,F)
        min_bounds, max_bounds = pmfutil.state_bounds_from_profile(xinterp,F)
        min_labels, max_labels = pmfutil.assign_state_labels(min_bounds,max_bounds)
        pmfutil.save_state_bounds(T,min_bounds,max_bounds,min_labels,max_labels)
        os.chdir("..")
        with open("%s_profile/%s_state_bounds.txt" % (coordname,T),"r") as fin:
            state_bins = np.array([ [float(x.split()[1]),float(x.split()[2])] 
                for x in fin.readlines() if x.split()[0].startswith("TS")])
    else:
        # Load state bounds
        with open("%s_profile/state_bounds.txt" % coordname,"r") as fin:
            state_bins = np.array([ [float(x.split()[1]),float(x.split()[2])] 
                for x in fin.readlines() if x.split()[0].startswith("TS")])

    # Calculate the contact probability in TS
    print "calculating TS"
    contact_args.bins = state_bins
    bin_edges, qi_vs_Q = simulation.calc.binned_contacts.calculate_binned_contacts_vs_q(contact_args)
    TS = qi_vs_Q[0,:]
    np.savetxt("binned_contacts_vs_%s/cont_prob_TS.dat" % coordname,TS)
    return TS

if __name__ == "__main__":
    variance = ["0.0001","0.0081","0.01","0.0625","0.49","0.64","0.81","1.00"]
    name = "PDZ"
    pairs = np.loadtxt("PDZca.contacts",dtype=int) - 1
    n_pairs = pairs.shape[0]
    n_residues = 95
    coordfile = "Qtanh_0_05.dat"
    coordname = coordfile.split(".")[0]
    tempsfile = "Tf_temps"
    tempsfile_backup = "ticatemps"

    parent_dirs = [ "random_native_%s" % variance[i] for i in range(len(variance)) ]
    sub_dirs = [ "replica_%d/%s/iteration_0" % (i,name) for i in range(1,11) ] 

    plot_replica_TS_grid(parent_dirs,sub_dirs,pairs,n_residues,coordfile,tempsfile,tempsfile_backup)


''' Plot a summary of MC2004 iteration


Description:


'''


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

from mpl_toolkits.axes_grid1 import make_axes_locatable

from bootFQ import get_F_with_error

colors = [('white')] + [(cm.Blues(i)) for i in xrange(1,256)]
global new_map
new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

colors2 = [('gray')] + [(cm.Blues(i)) for i in xrange(1,256)]
global new_map2
new_map2 = matplotlib.colors.LinearSegmentedColormap.from_list('new_map2', colors2, N=256)

#colors3 = [('gray')] + [(cm.Jet(i)) for i in xrange(1,256)]
#global new_map3
#new_map3 = matplotlib.colors.LinearSegmentedColormap.from_list('new_map3', colors3, N=256)

global SKIP_INTERACTIONS
SKIP_INTERACTIONS = [1,8,9]

def get_pairwise_params(pairwise_params_file,model_params_file):
    ''' Grab pairwise_params from file. '''
    model_param_values = np.loadtxt(model_params_file)
    model_params = []

    p_lines = [ x.rstrip("\n") for x in open(pairwise_params_file,"r").readlines() ]

    contacts = []
    pairwise_param_assignment = []
    pairwise_type = []
    pairwise_other_params = []

    for i in range(len(p_lines[1:])):
        data = p_lines[1+i].split()
        pairwise_param_assignment.append(int(data[2]))
        pairwise_type.append(int(data[3]))
        temp = []
        for otherparam in data[4:]:
            temp.append(float(otherparam))
        pairwise_other_params.append(temp)
        if not (int(data[3]) in SKIP_INTERACTIONS):
            contacts.append([int(data[0]),int(data[1])])
            if int(data[3]) in [3,5]:
                model_params.append(-model_param_values[i])
            else:
                model_params.append(model_param_values[i])

    contacts = np.array(contacts)
    pairwise_param_assignment = np.array(pairwise_param_assignment)
    pairwise_type = np.array(pairwise_type)
    model_params = np.array(model_params)

    return contacts,pairwise_param_assignment,model_params,pairwise_type,pairwise_other_params

def get_contact_probability(name,iteration,n_contacts,state_label,state_bound):

    cwd = os.getcwd()
    sub =  "%s/iteration_%d" % (name,iteration)
    os.chdir(sub)
    
    if not os.path.exists("cont_prob_%s.dat" % state_label):
        contact_probability = np.zeros(n_contacts,float) 
        n_frames = 0.
        temps = [ x.rstrip("\n") for x in open("long_temps_last", "r").readlines() ]
        for i in range(len(temps)):
            T = temps[i]
            Q = np.loadtxt("%s/Q.dat" % T)
            qimap = np.loadtxt("%s/qimap.dat" % T)

            state_indicator = ((Q > state_bound[0]).astype(int)*(Q < state_bound[1]).astype(int)).astype(bool)
            n_frames += float(sum(state_indicator.astype(int)))
            contact_probability += sum(qimap[(state_indicator == True),:])
        contact_probability /= n_frames
        np.savetxt("cont_prob_%s.dat" % state_label,contact_probability)
    else:
        contact_probability = np.loadtxt("cont_prob_%s.dat" % state_label)

    os.chdir(cwd)
    return contact_probability

def get_iteration_data(name,iteration):
    ''' Get summary data for iteration '''

    cwd = os.getcwd()
    sub =  "%s/iteration_%d" % (name,iteration)
    os.chdir(sub)

    Temps = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ]
    Tuse = open("long_temps_last","r").readlines()[0].rstrip("\n")
    if not os.path.exists("long_Tf"):
        Tf = float(Tuse.split("_")[0])
    else:
        Tf = float(open("long_Tf","r").read().rstrip("\n"))
    n_residues = len(open("%s/Native.pdb" % Temps[0],"r").readlines()) - 1

    contacts,pairwise_param_assignment,epsilons,pairwise_type,pairwise_other_params = get_pairwise_params("%s/pairwise_params" % Tuse,"%s/model_params" % Tuse)

    epsilon_map = np.zeros((n_residues,n_residues))
    for i in range(len(contacts)):
        epsilon_map[contacts[i,1]-1,contacts[i,0]-1] = epsilons[i] 
    
    Fs = []
    F_errs = []
    Bin_centers = []
    for i in range(len(Temps)): 
        filename = "%s/Q.dat" % Temps[i]
        F,F_err,bin_centers = get_F_with_error(filename,25,2,1)
        Fs.append(F)
        F_errs.append(F_err)
        Bin_centers.append(bin_centers)

    state_labels = []
    state_bounds = []
    for line in open("state_bounds.txt","r"):
        state_labels.append(line.split()[0])
        state_bounds.append([int(line.split()[1]),int(line.split()[2])])

    ddgsim = np.loadtxt("newton/sim_feature.dat")
    if not os.path.exists("newton/sim_feature_err.dat"):
        ddgsim_err = np.zeros(len(ddgsim))
    else:
        ddgsim_err = np.loadtxt("newton/sim_feature_err.dat")

    ddgexp = np.loadtxt("newton/target_feature.dat")
    if not os.path.exists("newton/target_feature_err.dat"):
        ddgexp_err = np.zeros(len(ddgsim))
    else:
        ddgexp_err = np.loadtxt("newton/target_feature_err.dat")

    ddGsim = np.zeros((len(ddgsim),2),float)
    ddGsim[:,0] = ddgsim
    ddGsim[:,1] = ddgsim_err
    ddGexp = np.zeros((len(ddgexp),2),float)
    ddGexp[:,0] = ddgexp
    ddGexp[:,1] = ddgexp_err

    n_contacts = len(contacts)

    os.chdir(cwd)

    return epsilons, epsilon_map, n_residues, contacts, n_contacts, Tf, state_labels, state_bounds, ddGsim, ddGexp, Fs, F_errs, Bin_centers

def plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues,individual=False):
    ''' Get  '''

    if os.path.exists("%s/epsilon_range" % name):
        temp = np.loadtxt("%s/epsilon_range" % name) 
        mineps = temp[0]
        maxeps = temp[1]
    else:
        mineps = min(epsilons) - 0.01
        maxeps = max(epsilons) + 0.01

    #plt.pcolor(epsilon_map,cmap=new_map2,vmin=mineps,vmax=maxeps)      ## Blues cmap
    plt.pcolor(epsilon_map,cmap="jet",vmin=mineps,vmax=maxeps)
    cbar = plt.colorbar()
    cbar.set_clim(mineps,maxeps)
    plt.xticks(range(0,n_residues+1,10))
    plt.yticks(range(0,n_residues+1,10))
    if individual:
        plt.title("%s params after iteration %d" % (name,iteration))
    plt.grid(True)
    plt.xlim(0,n_residues)
    plt.ylim(0,n_residues)

def plot_epsilon_histogram(name,iteration,epsilons,individual=False):
    plt.hist(epsilons,bins=20,histtype="stepfilled",facecolor=new_map(0.75),alpha=0.4)
    if individual:
        plt.title("%s iteration %d $\\overline{\\epsilon} = %.3f$    $\\sigma^2 = %.3f$ " % \
                (name,iteration,np.mean(epsilons),np.std(epsilons)))
    else:
        plt.title("$\\overline{\\epsilon} = %.3f$    $\\sigma^2 = %.3f$ " % \
                (np.mean(epsilons),np.std(epsilons)))
    plt.xlabel("$\\epsilon$")
    plt.ylabel("Frequency")
    plt.grid(True)

def plot_free_energy(name,iteration,n_contacts,Tf,Fs,F_errs,Bin_centers,state_labels,state_bounds,individual=False):

    colors = ["r","g","b","k","magenta","cyan"]
    barriers = []
    for i in range(len(Fs)):
        plt.plot(Bin_centers[i],Fs[i],color=colors[i],lw=2)
        plt.fill_between(Bin_centers[i],Fs[i] + F_errs[i],Fs[i] - F_errs[i],facecolor=colors[i],alpha=0.25)
        TS = ((Bin_centers[i] > state_bounds[1][0]).astype(int)*(Bin_centers[i] < state_bounds[1][1]).astype(int)).astype(bool)
        barr_height = max(Fs[i][TS])
        barriers.append(barr_height)
    barrier_height = max(barriers)
    plt.xlabel("Q")
    plt.ylabel("F(Q)")
    plt.xlim(0,n_contacts)
    plt.ylim(0,barrier_height+1)
    plt.vlines(state_bounds[0],0,1.5,lw=2)
    plt.vlines(state_bounds[1],barrier_height-1,barrier_height+0.5,lw=2)
    plt.vlines(state_bounds[2],0,1.5,lw=2)
    if individual:
        plt.title("%s iteration %d %.1f" % (name,iteration,Tf))
    else:
        plt.title("%.1f" % Tf)

def plot_ddG_comparison(name,iteration,ddGsim,ddGexp,individual=False):
    N = len(ddGsim)
    plt.errorbar(ddGsim[N/2:,0],ddGexp[N/2:,0],xerr=ddGsim[N/2:,1],yerr=ddGexp[N/2:,1],marker='o',linestyle="none",color='r',label="$\\Delta\\Delta G^{\\circ}$")
    plt.errorbar(ddGsim[:N/2,0],ddGexp[:N/2,0],xerr=ddGsim[:N/2,1],yerr=ddGexp[:N/2,1],marker='o',linestyle="none",color='b',label="$\\Delta\\Delta G^{\\dagger}$")
    plt.plot(range(-1,8),range(-1,8),color='k',lw=1)
    plt.xlabel("$\\Delta\\Delta G_{sim}$ (kT)")
    plt.ylabel("$\\Delta\\Delta G_{exp}$ (kT)")
    if individual:
        plt.title("%s iteration %d" % (name,iteration))
    maxddg = max([max(ddGsim[:,0]),max(ddGexp[:,0])])
    minddg = min([min(ddGsim[:,0]),min(ddGexp[:,0])])
    plt.xlim(minddg-0.5,maxddg+1)
    plt.ylim(minddg-0.5,maxddg+1)
    plt.grid(True)
    lg = plt.legend(loc=4)
    lg.draw_frame(False)

def plot_contact_probability(name,iteration,n_residues,contacts,state_label,state_bound,contact_probability,individual=False):
    ''' '''

    C = np.zeros((n_residues,n_residues),float)
    for j in range(len(contacts)):
        C[contacts[j,1]-1,contacts[j,0]-1] = contact_probability[j]

    #plt.pcolor(C,cmap=plt.get_cmap("Blues"),vmin=0.,vmax=1.0)
    plt.pcolor(C,cmap=new_map2,vmin=0.,vmax=1.0)
    if individual:
        cbar = plt.colorbar()
        cbar.set_clim(0,1)
    #for i in range(len(contacts)):
    #    plt.plot(contacts[i,0]-0.5,contacts[i,1]-0.5,marker='s',ms=6,markeredgecolor="k",markerfacecolor=new_map(contact_probability[i]))
    plt.xlim(0,n_residues)
    plt.ylim(0,n_residues)
    plt.xticks(range(0,n_residues,10))
    plt.yticks(range(0,n_residues,10))
    plt.grid(True)
    if individual:
        plt.title("%s iteration %d. %s contact probability" % (name,iteration,state_label))
    else:
        plt.title("%s contact probability" % (state_label))

def plot_contact_probability_subplot(name,iteration,n_residues,contacts,state_labels,Contact_maps):


    plt.figure(figsize=(9,4.2))
    for X in range(len(state_labels)):
        ax = plt.subplot(1,len(state_labels),X+1,aspect=1.0)
    
        C = np.zeros((n_residues,n_residues),float)
        for j in range(len(contacts)):
            C[contacts[j,1]-1,contacts[j,0]-1] = Contact_maps[X][j]

        #im = ax.pcolor(C,cmap=plt.get_cmap("Blues"),vmin=0.,vmax=1.0)
        im = ax.pcolor(C,cmap=new_map2,vmin=0.,vmax=1.0)
        plt.title("%s" % state_labels[X])
        #for i in range(len(contacts)):
        #    plt.plot(contacts[i,0]-0.5,contacts[i,1]-0.5,marker='s',ms=6,markeredgecolor="k",markerfacecolor=new_map(Contact_maps[X][i]))
        plt.xlim(0,n_residues)
        plt.ylim(0,n_residues)
        plt.xticks(range(0,n_residues,10))
        plt.yticks(range(0,n_residues,10))
        plt.grid(True)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_clim(0,1)

    plt.suptitle("%s iteration %d" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/contact_prob_all.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/contact_prob_all.png" % (name,iteration))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration

    epsilons, epsilon_map, n_residues, contacts, n_contacts, Tf, state_labels, state_bounds, ddGsim, ddGexp, Fs, F_errs, Bin_centers = get_iteration_data(name,iteration)
    
    print "Plotting summary for %s iteration %d..." % (name,iteration)

    if not os.path.exists("%s/iteration_%d/summary" % (name,iteration)):
        os.mkdir("%s/iteration_%d/summary" % (name,iteration))

    Contact_maps = []
    for X in range(len(state_labels)):
        plt.figure()
        contact_probability = get_contact_probability(name,iteration,n_contacts,state_labels[X],state_bounds[X])
        plot_contact_probability(name,iteration,n_residues,contacts,state_labels[X],state_bounds[X],contact_probability,individual=True)
        Contact_maps.append(contact_probability)
        print " Saving: %s/iteration_%d/contact_prob_%s.pdf          - %s contact probabilities" % (name,iteration,state_labels[X],state_labels[X])
        plt.savefig("%s/iteration_%d/summary/contact_prob_%s.pdf" % (name,iteration,state_labels[X]))
        plt.savefig("%s/iteration_%d/summary/contact_prob_%s.png" % (name,iteration,state_labels[X]))
        plt.close()

    print " Saving subplot: %s/iteration_%d/summary/contact_prob_all.pdf       - contact probabilities" % (name,iteration)
    plot_contact_probability_subplot(name,iteration,n_residues,contacts,state_labels,Contact_maps)
    
    print "  Saving: %s/iteration_%d/summary/current_epsilon_map.pdf    - epsilon map" % (name,iteration)
    plt.figure()
    plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues,individual=True)
    plt.savefig("%s/iteration_%d/summary/current_epsilon_map.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/current_epsilon_map.png" % (name,iteration))
    plt.close()

    print "  Saving: %s/iteration_%d/summary/current_epsilon_hist.pdf   - epsilon histogram" % (name,iteration)
    plt.figure()
    plot_epsilon_histogram(name,iteration,epsilons,individual=True)
    plt.savefig("%s/iteration_%d/summary/current_epsilon_hist.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/current_epsilon_hist.png" % (name,iteration))
    plt.close()

    print "  Saving: %s/iteration_%d/summary/FreeEnergy_Q.pdf            - free energy" % (name,iteration)
    plt.figure()
    plot_free_energy(name,iteration,n_contacts,Tf,Fs,F_errs,Bin_centers,state_labels,state_bounds,individual=True)
    plt.savefig("%s/iteration_%d/summary/FreeEnergy_Q.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/FreeEnergy_Q.png" % (name,iteration))
    plt.close()

    print "  Saving: %s/iteration_%d/compareddG.pdf              - ddG comparison" % (name,iteration)
    plt.figure()
    plot_ddG_comparison(name,iteration,ddGsim,ddGexp,individual=True)
    plt.savefig("%s/iteration_%d/summary/compareddG.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/summary_%s_%d.png" % (name,iteration,name,iteration))
    plt.close()

    print " Summary figure: %s/iteration_%d/summary/summary_%s_%d.pdf" % (name,iteration,name,iteration)
    plt.figure(figsize=(12,10))
    plt.subplot(2,2,1)
    plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues)
    plt.subplot(2,2,2)
    plot_epsilon_histogram(name,iteration,epsilons)
    plt.subplot(2,2,3)
    plot_free_energy(name,iteration,n_contacts,Tf,Fs,F_errs,Bin_centers,state_labels,state_bounds)
    plt.subplot(2,2,4)
    plot_ddG_comparison(name,iteration,ddGsim,ddGexp)
    plt.suptitle("%s iteration %d" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/summary_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/summary/summary_%s_%d.png" % (name,iteration,name,iteration))
    #plt.show()

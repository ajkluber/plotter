import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

from plotter.cube_cmap import cubecmap

class Dataset(object):
    """Plotting one thing against a reaction coordinate"""
    def __init__(self, parent_dirs, sub_dirs, datafile, binfile, coordname, 
            xlabels, ylabels, titles, saveas, save_formats, xytext=(0.1, 0.1),
            xscale_factor=1., yscale_factor=1.):
        self.parent_dirs = parent_dirs
        self.sub_dirs = sub_dirs
        self.datafile = datafile
        self.binfile = binfile
        self.coordname = coordname
        self.xlabels = xlabels
        self.ylabels = ylabels
        self.titles = titles
        self.xytext = xytext 
        self.xscale_factor = xscale_factor
        self.yscale_factor = yscale_factor
        self.saveas = saveas
        self.save_formats = save_formats
        self.rep_savedir = os.path.dirname(self.datafile)
        self.rep_avg_saveas = os.path.basename(self.datafile)
        self.rep_std_saveas = self.rep_avg_saveas.split(".dat")[0] + "_std.dat"
        self.dimension = 1

    def get_data(self):
        """Get data from subdirectory"""
        data = self.yscale_factor*np.loadtxt(self.datafile, dtype=float)
        mid_bin = self.xscale_factor*np.loadtxt(self.binfile, dtype=float)
        if len(data.shape) == 2:
            self.dimension = data.shape[1]
        else:
            self.dimension = 1
        return mid_bin, data

    def save_rep_avg(self, interp_points, rep_bin_avg, rep_bin_std, rep_min_max):
        """Save replica averaged data"""
        cwd = os.getcwd()
        if not os.path.exists(self.rep_savedir):
            os.mkdir(self.rep_savedir)
        os.chdir(self.rep_savedir)
        np.savetxt(self.rep_avg_saveas, rep_bin_avg/self.yscale_factor)
        np.savetxt(self.rep_std_saveas, rep_bin_avg/self.yscale_factor)
        np.savetxt("mid_bin.dat", interp_points/self.xscale_factor)
        np.savetxt("states.dat", rep_min_max/self.xscale_factor)
        os.chdir(cwd)
    
def plot_replica_grid_with_colors(dataset, xlabel, ylabel, title, annotations,
            coloridxs, grid_dims, saveas, saveas_formats=["png","pdf"], display=False, y_ticks=None):
    """ Plot replica observables with their replica average  NOT DONE"""
    cwd = os.getcwd()

    parent_dirs = dataset.parent_dirs
    sub_dirs = dataset.sub_dirs
    n_dirs = len(parent_dirs)
     
    all_interp_points, all_rep_avg_data, all_min_max = get_replica_average(dataset)
    fig, axes = plt.subplots(grid_dims[0], grid_dims[1], sharex=True, sharey=True)
    counter = 0
    for i in range(n_dirs):
         Plot individual replicas
        rep_bins, rep_data = get_replicas(dataset, i)
        ax = axes[counter / grid_dims[1], counter % grid_dims[1]]
        for rep_idx in range(len(rep_bins)):
            ax.plot(rep_bins[rep_idx], rep_data[rep_idx], color=cubecmap(coloridxs[i]), alpha=0.3)

        # Plot bold replica average
        ax.plot(all_interp_points[i], all_rep_avg_data[i], lw=3, color=cubecmap(coloridxs[i]))
        counter += 1

        if y_ticks is not None:
            ax.set_yticks(y_ticks)
            ax.set_ylim(y_ticks[0],y_ticks[-1])

        ax.annotate("%s" % annotations[i], xy=(0.5, 0.5), xytext=dataset.xytext,
                            xycoords="axes fraction", textcoords="axes fraction")

    # Rotate ticklabels.
    for j in range(axes.shape[1]): 
        plt.setp(axes[-1,j].xaxis.get_majorticklabels(), rotation=270)

    big_ax = fig.add_subplot(111)
    big_ax.grid(False)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.set_ylabel(ylabel)
    big_ax.set_xlabel(xlabel)
    fig.suptitle(title, fontsize=18)
    fig.subplots_adjust(hspace=0, wspace=0)

    os.chdir("%s" % cwd)
    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    for format in saveas_formats:
        fig.savefig("%s.%s" % (saveas, format), bbox_inches="tight")
    if display:
        plt.show()
    plt.close(fig)
    os.chdir(cwd)

def plot_replica_grid(dataset, display=False):
    """ Plot grid of some observable for replicas"""
    cwd = os.getcwd()
    all_rep_avg = []
    parent_dirs = dataset.parent_dirs
    sub_dirs = dataset.sub_dirs
    for i in range(len(parent_dirs)):
        # Plot all replicas for a given parent directory
        print parent_dirs[i]
        fig1,axes = plt.subplots(3,4,sharex=True,sharey=True)
        min_rep = []
        max_rep = []
        counter = 0
        for rep_idx in range(len(sub_dirs)):
            # Plot dataset
            os.chdir("%s/%s/%s" % (cwd,parent_dirs[i],sub_dirs[rep_idx]))
            ax = axes[counter / 4, counter % 4]

            if hasattr(dataset,"plot"):
                dataset.plot(ax)
            else:
                mid_bin, data = dataset.get_data()
                ax.plot(mid_bin, data)

            ax.annotate("replica %d" % (rep_idx + 1), xy=(0.5, 0.5), xytext=dataset.xytext, 
                                xycoords="axes fraction", textcoords="axes fraction")

            # Get state locations
            min_rep.append(dataset.xscale_factor*np.loadtxt("%s_profile/minima.dat" % dataset.coordname))
            max_rep.append(dataset.xscale_factor*np.loadtxt("%s_profile/maxima.dat" % dataset.coordname))
            counter += 1

        # Rotate ticklabels.
        for j in range(axes.shape[1]): 
            plt.setp(axes[-1,j].xaxis.get_majorticklabels(), rotation=270)

        # Fill regions for three states.
        counter = 0
        for rep_idx in range(len(sub_dirs)):
            ax = axes[counter / 4, counter % 4]
            fill_areas = [ np.array([p - 4, p + 4]) for p in min_rep[rep_idx] ]
            fill_areas += [ np.array([max_rep[rep_idx] - 4, max_rep[rep_idx] + 4]) ]
            y0,y1 = ax.get_ylim()
            for m in range(len(fill_areas)):
                ax.fill_between(fill_areas[m],y0,y1,alpha=0.5,color='#5DA5DA')
            counter += 1

        # common xlabel/ylabel and title
        big_ax = fig1.add_subplot(111)
        big_ax.set_axis_bgcolor('none')
        big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        big_ax.set_ylabel(dataset.ylabels[i])
        big_ax.set_xlabel(dataset.xlabels[i])
        fig1.suptitle(dataset.titles[i],fontsize=18)
        fig1.subplots_adjust(hspace=0,wspace=0)

        # save figure
        os.chdir("%s/%s" % (cwd,parent_dirs[i]))
        if not os.path.exists("plots"):
            os.mkdir("plots")
        os.chdir("plots")
        for format in dataset.save_formats:
            fig1.savefig("%s.%s" % (dataset.saveas,format),bbox_inches="tight")
        if display:
            plt.show()
        plt.close(fig1)
        os.chdir(cwd)

def plot_replica_average(dataset, xlabel, ylabel, title, annotations, coloridx, display=False):
    """ Plot replica averaged observable"""
    cwd = os.getcwd()
    ninterp = 500
    parent_dirs = dataset.parent_dirs
    sub_dirs = dataset.sub_dirs
    sqr_num = int(np.ceil(np.sqrt(float(len(parent_dirs)))))
    fig1, axes = plt.subplots(sqr_num, sqr_num, sharex=True, sharey=True)
    counter = 0
    # Plot replica averaged quantity 
    all_interp_points, all_rep_avg_data, all_min_max = get_replica_average(dataset)

    # Plot replica averaged curves in a grid.
    for i in range(len(parent_dirs)):
        interp_points = all_interp_points[i]
        rep_avg_data = all_rep_avg_data[i]
        ax = axes[counter / sqr_num, counter % sqr_num]
        if hasattr(dataset, "plot_rep_avg"):
            dataset.plot_rep_avg(ax, interp_points, rep_avg_data)
        else:
            ax.plot(interp_points, rep_avg_data, color=cubecmap(coloridx[i]))

        if hasattr(dataset,"ylim"):
            ax.set_ylim(*dataset.ylim)

        ax.annotate("%s" % annotations[i], xy=(0.5, 0.5), xytext=dataset.xytext,
                            xycoords="axes fraction", textcoords="axes fraction")
        counter += 1

    # Plot average position of basins.
    counter = 0
    for i in range(len(parent_dirs)):
        y0,y1 = ax.get_ylim()
        ax = axes[counter / sqr_num, counter % sqr_num]
        ax.fill_between(np.array([all_min_max[i][0] - 5, all_min_max[i][0] + 5]),y0,y1,alpha=0.5,color='#5DA5DA')
        ax.fill_between(np.array([all_min_max[i][1] - 5, all_min_max[i][1] + 5]),y0,y1,alpha=0.5,color='#5DA5DA')
        ax.fill_between(np.array([all_min_max[i][1] - 5, all_min_max[i][2] + 5]),y0,y1,alpha=0.5,color='#5DA5DA')
        counter += 1

    # Rotate ticklabels.
    for j in range(axes.shape[1]): 
        plt.setp(axes[-1,j].xaxis.get_majorticklabels(), rotation=270)

    # common xlabel/ylabel and title
    big_ax = fig1.add_subplot(111)
    big_ax.grid(False)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.set_ylabel(ylabel)
    big_ax.set_xlabel(xlabel)
    big_ax.set_title(title)
    fig1.subplots_adjust(hspace=0,wspace=0)

    # plot together on one axis as well.
    if dataset.dimension == 1:
        fig2,ax2 = plt.subplots(1,1)
        for i in range(len(parent_dirs)):
            avg_data = all_rep_avg_data[i]
            interp_points = all_interp_points[i]
            ax2.plot(interp_points, avg_data, color=cubecmap(coloridx[i]),label=annotations[i])

            tempU = (interp_points - all_min_max[i][0])**2
            U = np.where(tempU == tempU.min())[0]
            ax2.plot(interp_points[U], avg_data[U], marker='o', ms=8, color=cubecmap(coloridx[i]))

            tempTS = (interp_points - all_min_max[i][1])**2
            TS = np.where(tempTS == tempTS.min())[0]
            ax2.plot(interp_points[TS], avg_data[TS], marker='s', ms=8, color=cubecmap(coloridx[i]))

            tempN = (interp_points - all_min_max[i][2])**2
            N = np.where(tempN == tempN.min())[0]
            ax2.plot(interp_points[N], avg_data[N], marker='>', ms=8, color=cubecmap(coloridx[i]))

            if hasattr(dataset,"ylim"):
                ax2.set_ylim(*dataset.ylim)
        
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax2.set_ylabel(ylabel)
        ax2.set_xlabel(xlabel)
        ax2.set_title(title)

    # save figure
    os.chdir(cwd)
    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    for format in dataset.save_formats:
        fig1.savefig("%s.%s" % (dataset.saveas,format),bbox_inches="tight")
        if dataset.dimension == 1:
            fig2.savefig("one_%s.%s" % (dataset.saveas,format),bbox_inches="tight")
    if display:
        plt.show()
    plt.close(fig1)


def plot_protein_grid(names, datasets, xlabel, ylabel, title, coloridxs, 
         saveas, save_formats=["png","pdf"], display=False, x_ticks=None, y_ticks=None, aspect_equal=None):
    """Plot all proteins on one grid"""

    fig, axes = plt.subplots(1, len(names), sharey=True, figsize=(13,4.))
    

    cwd = os.getcwd()
    state_markers = ["o","s",">"]
    for n in range(len(names)):
        dataset = datasets[n]
        coloridx = coloridxs[n]
        name = names[n]
        os.chdir(name)
        all_interp_points, all_rep_avg_data, all_min_max = get_replica_average(dataset)
        os.chdir("..")

        ax = axes[n]
        # Loop over variances
        for i in range(len(all_rep_avg_data)):
            avg_data = all_rep_avg_data[i]
            interp_points = all_interp_points[i]
            ax.plot(interp_points, avg_data, color=cubecmap(coloridx[i]))

            # Indicate location of minima, maxima 
            for s in range(len(all_min_max[i])):
                temp = (interp_points - all_min_max[i][s])**2
                state = np.where(temp == temp.min())[0]
                ax.plot(interp_points[state], avg_data[state], 
                    marker=state_markers[s], ms=8, color=cubecmap(coloridx[i]))

            if hasattr(dataset,"ylim"):
                ax.set_ylim(*dataset.ylim)

        ax.annotate(name, xy=(0.5, 0.5), xytext=(0.05,0.92), fontsize=20,
                            xycoords="axes fraction", textcoords="axes fraction")
        if aspect_equal:
            ax.set_aspect('equal', 'datalim')

        if y_ticks is not None:
            ax.set_yticks(y_ticks)
            ax.set_ylim(y_ticks[0],y_ticks[-1])
        if x_ticks is not None:
            ax.set_xticks(x_ticks)
            ax.set_xlim(x_ticks[0],x_ticks[-1])

    # Rotate ticklabels.
    for j in range(len(axes)): 
        plt.setp(axes[j].xaxis.get_majorticklabels(), rotation=270)

    # common xlabel/ylabel and title
    big_ax = fig.add_subplot(111)
    big_ax.grid(False)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    #big_ax.set_ylabel(ylabel)
    big_ax.text(-.085, 0.9, ylabel, rotation="vertical", fontsize=16, transform=big_ax.transAxes)
    big_ax.set_xlabel(xlabel)
    fig.suptitle(title,fontsize=20)
    fig.subplots_adjust(hspace=0,wspace=0)

    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    # save figures
    for format in save_formats:
        fig.savefig("%s.%s" % (saveas,format), bbox_inches="tight")
    if display:
        plt.show() 
    os.chdir("..")

def get_replicas(dataset, i):
    cwd = os.getcwd()
    parent_dirs = dataset.parent_dirs
    sub_dirs = dataset.sub_dirs

    rep_data = []
    rep_bins = []
    # Get replica data
    for rep_idx in range(len(sub_dirs)):
        os.chdir("%s/%s/%s" % (cwd, parent_dirs[i], sub_dirs[rep_idx]))
        mid_bin, data = dataset.get_data()
        rep_data.append(data)
        rep_bins.append(mid_bin)
    os.chdir(cwd)

    return rep_bins, rep_data

def get_replica_average(dataset):
    """Collect replica averaged observable"""
    cwd = os.getcwd()
    ninterp = 500
    parent_dirs = dataset.parent_dirs
    sub_dirs = dataset.sub_dirs
    all_min_max = []
    all_interp_points = []
    all_rep_avg_data = []
    for i in range(len(parent_dirs)):
        print parent_dirs[i]
        if os.path.exists("%s/%s/%s/%s" % (cwd, parent_dirs[i], dataset.rep_savedir, dataset.rep_avg_saveas)):
            # Load data if it has been calculated.
            os.chdir("%s/%s/%s" % (cwd, parent_dirs[i], dataset.rep_savedir))
            rep_avg_data = dataset.yscale_factor*np.loadtxt(dataset.rep_avg_saveas)
            rep_std_data = dataset.yscale_factor*np.loadtxt(dataset.rep_std_saveas)
            interp_points = dataset.xscale_factor*np.loadtxt("mid_bin.dat")
            rep_min_max = dataset.xscale_factor*np.loadtxt("states.dat")
            os.chdir(cwd)
        else:
            rep_data = []
            rep_bins = []
            rep_interp = []
            min_rep = []
            max_rep = []
            # Calculate replica average.
            for rep_idx in range(len(sub_dirs)):
                os.chdir("%s/%s/%s" % (cwd, parent_dirs[i], sub_dirs[rep_idx]))
                mid_bin, data = dataset.get_data()

                #print "", mid_bin
                if dataset.dimension > 1:
                    interp = []
                    for n in range(dataset.dimension):
                        interp.append(interp1d(mid_bin, data[:,n]))
                else:
                    interp = interp1d(mid_bin, data, kind="cubic")
                rep_data.append(data)
                rep_bins.append(mid_bin)
                rep_interp.append(interp)

                # Get state locations
                min_rep.append(dataset.xscale_factor*np.loadtxt("%s_profile/minima.dat" % dataset.coordname, dtype=float))
                max_rep.append(dataset.xscale_factor*np.loadtxt("%s_profile/maxima.dat" % dataset.coordname, dtype=float))

            # Average positions of the minima and maxima.
            avgU = np.mean([ min_rep[r][0] for r in range(len(min_rep)) ])
            avgN = np.mean([ min_rep[r][1] for r in range(len(min_rep)) ])
            avgTS = np.mean([ max_rep[r] for r in range(len(max_rep)) ])
            rep_min_max = np.array([avgU, avgTS, avgN])

            # Find common domain to interpolate data over.
            common_min = np.max([ mid_bin.min() for mid_bin in rep_bins ])
            common_max = np.min([ mid_bin.max() for mid_bin in rep_bins ])
            interp_points = np.linspace(common_min, common_max, ninterp)

            # Interpolate data
            if dataset.dimension > 1:
                new_rep_data = np.zeros((len(rep_interp), dataset.dimension, ninterp), float)
                for r in range(len(rep_interp)):
                    for n in range(dataset.dimension):
                        new_rep_data[r,n,:] = rep_interp[r][n](interp_points)
            else:
                new_rep_data = np.array([ rep_interp[r](interp_points) for r in range(len(rep_interp)) ])

            # Calculate replica average
            rep_avg_data = np.mean(new_rep_data, axis=0)
            rep_std_data = np.std(new_rep_data, axis=0)

            # Save replica-average
            os.chdir("%s/%s" % (cwd, parent_dirs[i]))
            dataset.save_rep_avg(interp_points, rep_avg_data, rep_std_data, rep_min_max)

        all_interp_points.append(interp_points)
        all_rep_avg_data.append(rep_avg_data)
        all_min_max.append(rep_min_max)
    os.chdir(cwd)

    return all_interp_points, all_rep_avg_data, all_min_max

if __name__ == "__main__":
    pass

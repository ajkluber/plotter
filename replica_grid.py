import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

#import plotter.pmfutil as pmfutil
#from plotter.cube_cmap import cubecmap

class Dataset(object):
    """Plotting one thing against a reaction coordinate"""
    def __init__(self, parent_dirs, sub_dirs, datafile, binfile, coordname, 
            xlabels, ylabels, titles, saveas, save_formats, xytext=(0.1, 0.1)):
        self.parent_dirs = parent_dirs
        self.sub_dirs = sub_dirs
        self.datafile = datafile
        self.binfile = binfile
        self.coordname = coordname
        self.xlabels = xlabels
        self.ylabels = ylabels
        self.titles = titles
        self.xytext = xytext 
        self.saveas = saveas
        self.save_formats = save_formats

    def get_data(self):
        data = np.loadtxt(self.datafile)
        mid_bin = np.loadtxt(self.binfile)
        return mid_bin, data
    
def plot_replica_grid(dataset, display=False):
    """ Plot grid of some observable contact probabilities for replicas

    """
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
            mid_bin, data = dataset.get_data()
            ax.plot(mid_bin, data)
            ax.annotate("replica %d" % (rep_idx + 1), xy=(0.5, 0.5), xytext=dataset.xytext, 
                                xycoords="axes fraction", textcoords="axes fraction")

            # Get state locations
            min_rep.append(np.loadtxt("%s_profile/minima.dat" % dataset.coordname))
            max_rep.append(np.loadtxt("%s_profile/maxima.dat" % dataset.coordname))
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

if __name__ == "__main__":
    pass

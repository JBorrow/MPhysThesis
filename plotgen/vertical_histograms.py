""" Makes the vertical histograms thing on the cluster """

import survis
import os
import pickle
import numpy as np
from tqdm import tqdm
from scipy.optimize import curve_fit
import matplotlib

matplotlib.use('TkAgg')

# ---- Global Parameters
snapshot = 97
bbox_x = [-40, 40]
bbox_y = bbox_x
res_elem_lr = 0.5
res_elem_nr = 0.2
res_lr = survis.helper.get_res(res_elem_lr, bbox_x, bbox_y)
res_nr = survis.helper.get_res(res_elem_nr, bbox_x, bbox_y)

names = ['custom', 'custom_highres', 'custom_lowres']
# ----------------------


class Vertical(object):
    def __init__(self, filename, res, bbox_x, bbox_y, max_rad=30, n_bins=5, sf=1):

        data = survis.preprocess.DataGridder(filename,
                                             res[0],
                                             res[1],
                                             bbox_x[0],
                                             bbox_x[1],
                                             bbox_y[0],
                                             bbox_y[1],
                                             False)

        self.filename = filename
        self.res = res
        self.bbox_x = bbox_x
        self.bbox_y = bbox_y
        self.max_rad = max_rad
        self.n_bins = n_bins
        self.sf = sf  #scale factor

        self.z_data = data.gas['Coordinates'][:, 2]

        self.popt, self.pcov = self.get_scale_height(self.z_data, nbins=1000)

        return


    def bin_cent(self, bins):
        return 0.5*(bins[1:] + bins[:-1])


    def to_fit(self, z, norm, Z, offset):
        return norm*(1/(np.cosh((z+offset)/Z)**2))


    def get_scale_height(self, z_data, min=-0.2, max=0.2, nbins=25):
        # First we must bin the data

        def to_fit(z, norm, Z, offset):
            return norm*(1/(np.cosh((z+offset)/Z)**2))

        bins = np.arange(min, max, self.sf*(max-min)/nbins)
        n, bins = np.histogram(z_data, bins)
        bincenters = self.bin_cent(bins)

        print(bincenters)

        return curve_fit(to_fit, bincenters, n)


    def make_histogram(self, ax, nbins=1000, max=0.2, min=-0.2, name="", namemax=800):
        ax.hist(self.z_data, bins=int(1000/self.sf), range=(min, max))

        x_data = np.arange(min, max, (max - min)/nbins)
        y_data = self.to_fit(x_data, self.popt[0], self.popt[1], self.popt[2])
        ax.plot(x_data, y_data)

        ax.text(min+0.02, namemax-100, "$Z = {:2.1f}$ pc\n{}".format(1000*self.popt[1], name), fontsize=8)

        ax.set_xlim(min, max)

        return ax



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    # Run in scripting mode, and actually perform the analysis
    output = {}

    for simulation in tqdm(names):
        filename = "{}/snapshot_{:03d}.hdf5".format(simulation, snapshot)
        tqdm.write("Attempting to analyse {} {}".format(simulation, snapshot))

        if "lowres" in simulation:
            data = Vertical(filename, res_lr, bbox_x, bbox_y, sf=10)
        else:
            data = Vertical(filename, res_nr, bbox_x, bbox_y)
        tqdm.write("Analysis of {} {} successful.".format(simulation, snapshot))
    
        output[simulation] = data

    
    # Now we do the postprocessing

    f, axes = plt.subplots(1, 3, sharey=True, figsize=(6.3, 3)) 

    for ax, simulation in zip(axes, output):
        data = output[simulation]

        # Issue; plot overlaps with labels
        if simulation != "custom":
            label = simulation[7:]
        else:
            label = simulation

        data.make_histogram(ax, name=label)

        ax.set_ylim(0, 800)


    # Remove extra stuff

    nbins = len(axes[0].get_xticklabels())
    axes[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
    axes[1].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
    axes[2].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))

    axes[1].set_xlabel("Height above midplane [kpc]")
    axes[0].set_ylabel("Number of particles in bin")

    plt.tight_layout()
    f.subplots_adjust(wspace=0)

    f.savefig('v_height_hist.pdf', dpi=600)

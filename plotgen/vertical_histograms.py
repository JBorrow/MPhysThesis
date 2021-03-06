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
    def __init__(self, filename, res, bbox_x, bbox_y, max_rad=30, n_bins=5, sf=1, exclude=0):

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

        rads = np.sqrt(np.sum(np.square(data.gas['Coordinates']), 1))
        rads_star = np.sqrt(np.sum(np.square(data.star['Coordinates']), 1))

        mask = (rads < exclude)
        mask_star = (rads_star < exclude)

        z_data = data.gas['Coordinates'][:, 2]
        z_data_stars = data.star['Coordinates'][:, 2]

        mx = np.ma.masked_array(z_data, mask)
        mx_star = np.ma.masked_array(z_data_stars, mask_star)

        self.z_data = mx.compressed()
        self.z_data_stars = mx_star.compressed()

        self.popt, self.pcov = self.get_scale_height(self.z_data, nbins=1000)
        self.popt_star, self.pcov_star = self.get_scale_height(self.z_data_stars, nbins=1000)

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
        x_data = np.arange(min, max, (max - min)/nbins)

        for z, popt in zip([self.z_data, self.z_data_stars], [self.popt, self.popt_star]):
            ax.hist(z, bins=int(1000/self.sf), range=(min, max), alpha=0.5, histtype='stepfilled')
            y_data = self.to_fit(x_data, popt[0], popt[1], popt[2])
            ax.plot(x_data, y_data)

        ax.text(min+0.02, namemax, "$Z_g = {:2.1f}$ pc\n$Z_* = {:3.0f}$ pc\n{}".format(1000*abs(self.popt[1]), abs(self.popt_star[1]*1000), name), fontsize=8)

        ax.set_xlim(min, max)

        return ax



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    # Run in scripting mode, and actually perform the analysis
    output = {}
    output_excl = {}  # Has the particles that are more than 0.5kpc out of the center.

    for simulation in tqdm(names):
        filename = "{}/snapshot_{:03d}.hdf5".format(simulation, snapshot)
        tqdm.write("Attempting to analyse {} {}".format(simulation, snapshot))

        if "lowres" in simulation:
            data = Vertical(filename, res_lr, bbox_x, bbox_y, sf=10)
        else:
            data = Vertical(filename, res_nr, bbox_x, bbox_y)
        tqdm.write("Analysis of {} {} successful.".format(simulation, snapshot))
    
        output[simulation] = data


    for simulation in tqdm(names):
        filename = "{}/snapshot_{:03d}.hdf5".format(simulation, snapshot)
        tqdm.write("Attempting to analyse {} {}".format(simulation, snapshot))

        if "lowres" in simulation:
            data = Vertical(filename, res_lr, bbox_x, bbox_y, sf=10, exclude=0.5)
        else:
            data = Vertical(filename, res_nr, bbox_x, bbox_y, exclude=0.5)
        tqdm.write("Analysis of {} {} successful.".format(simulation, snapshot))
    
        output_excl[simulation] = data

    
    # Now we do the postprocessing

    f, axes = plt.subplots(2, 3, figsize=(6.3, 6)) 

    for ax, simulation in zip(axes[0, :], output):
        data = output[simulation]

        # Issue; plot overlaps with labels
        if simulation != "custom":
            label = simulation[7:]
        else:
            label = simulation

        data.make_histogram(ax, name=label, namemax=660)

        ax.set_ylim(0, 800)


    for ax, simulation in zip(axes[1, :], output_excl):
        data = output_excl[simulation]

        # Issue; plot overlaps with labels
        if simulation != "custom":
            label = simulation[7:]
        else:
            label = simulation

        data.make_histogram(ax, name=label, namemax=82)

        ax.set_ylim(0, 100)


    axes[0, 0].tick_params(labelbottom='off')
    axes[0, 1].tick_params(labelbottom='off', labelleft='off')
    axes[0, 2].tick_params(labelbottom='off', labelleft='off')
    axes[1, 1].tick_params(labelleft='off')
    axes[1, 2].tick_params(labelleft='off')

    # Remove extra stuff

    nbins = len(axes[1, 0].get_xticklabels())
    nbinsy = len(axes[1, 0].get_yticklabels())
    axes[1, 0].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
    axes[1, 1].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
    axes[1, 2].xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
    axes[1, 0].yaxis.set_major_locator(MaxNLocator(nbins=nbinsy, prune='upper'))

    axes[1, 1].set_xlabel("Height above midplane [kpc]")
    axes[0, 0].set_ylabel("Number of particles in bin")

    plt.tight_layout()
    f.subplots_adjust(wspace=0, hspace=0)

    f.savefig('v_height_hist.pdf', dpi=600)

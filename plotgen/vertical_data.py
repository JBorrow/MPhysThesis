""" Gets the data for the vertical table thing on the cluster """

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

names = ['custom', 'custom_highres', 'custom_lowres', 'default', 'default_lowres']

G = 4.302e-3
# ----------------------

def disp(sg, fg=0.1):
    """ The dispersion from the Martizzi model in km/s """
    pfm = 300000
    f = 0.4
    F = 0.5

    return 1.8 * (f/F)**(3/5) * G**(2/5) * pfm**(1/5) * fg**(-2/5) * sg**(1/5)


def tau_sfr(rhog, fg=0.1):
    f = 0.4
    F = 0.5
    pfm = 3000
    pf = pfm*100

    prefactor = 0.45 * (fg**(3./2) * F**(7/4))/(np.pi * f**(3/4) * G**(7/8) * pf**(1/4))

    conversion = 3.171e-8 * 3.086e13

    return prefactor * pfm * rhog**(-5/8) * conversion


class Vertical(object):
    def __init__(self, filename, res, bbox_x, bbox_y, max_rad=30, n_bins=5, sf=1, inner=0, outer=0.5):

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
        self.inner = inner
        self.outer = outer

        rads = np.sqrt(np.sum(np.square(data.gas['Coordinates']), 1))
        rads_star = np.sqrt(np.sum(np.square(data.star['Coordinates']), 1))

        mask = np.logical_or((rads < inner), (rads > outer))
        mask_star = np.logical_or((rads_star < inner), (rads_star > outer))

        dens = data.gas['Density']
        self.gas_mass_res = data.gas_mass

        z_data = data.gas['Coordinates'][:, 2]
        z_data_stars = data.star['Coordinates'][:, 2]

        mx = np.ma.masked_array(z_data, mask)
        mx_star = np.ma.masked_array(z_data_stars, mask_star)
        m_dens = np.ma.masked_array(dens, mask)

        self.z_data = mx.compressed()
        self.z_data_stars = mx_star.compressed()
        self.dens = m_dens.compressed()

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
    

    def local_density_mean(self):
        """ Gets the local mean density in nh cm3 """

        # the simulation units are msun / kpc ^3
        local = np.mean(self.dens)

        return local


    def expected_lj_mart(self, sg, rhog):
        disp_val = disp(sg)

        lj =  disp_val * np.sqrt(9/(8*G)) / np.sqrt(rhog)

        return lj


    def local_surface_density(self):
        n = len(self.z_data)

        area = np.pi * 2 * (self.outer**2 - self.inner**2) * 1e6

        return n * self.gas_mass_res / area


    def print_data(self, latex=False, name=None):
        """ Prints the data for making the table """
        factor_3d_dens = 4.04367e-8

        rho_g = self.local_density_mean()
        sg = self.local_surface_density()

        expected_lj = (sg/(rho_g))*1e9  # in pc
        proj_lj = self.expected_lj_mart(sg, rho_g/1e9)

        if not latex:
            print("rho: {} nh/cm3, sg: {} ms/pc2, sg/rho: {} pc, Z: {} pc, lj: {} pc, sfr: {} gyr".format(
                rho_g * factor_3d_dens,
                sg,
                expected_lj,
                abs(self.popt[1]) * 1e3,
                proj_lj,
                tau_sfr(rho_g/1e9)/1e9))
        else:
            print("{} & {:.3g} & {:.3g} & {:.3g} & {:.3g} &  {:.3g} \\\\".format(
                name,
                rho_g * factor_3d_dens,
                sg,
                expected_lj,
                proj_lj,
                abs(self.popt[1])*1e3))


        return



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import sys
    # Run in scripting mode, and actually perform the analysis
    output = {}
    output_excl = {}  # Has the particles that are more than 0.5kpc out of the center.

    latex = "--latex" in sys.argv


    for simulation in names:
        filename = "{}/snapshot_{:03d}.hdf5".format(simulation, snapshot)
        if not latex:
            tqdm.write("Attempting to analyse {} {}".format(simulation, snapshot))

        if "lowres" in simulation:
            data = Vertical(filename, res_lr, bbox_x, bbox_y, sf=10, inner=0, outer=0.5)
        else:
            data = Vertical(filename, res_nr, bbox_x, bbox_y, inner=0, outer=0.5)

        name = r"{{\tt {}}} & $r<0.5$".format(simulation)
        data.print_data(latex, name)
    
        output[simulation] = data

    
    if not latex:
        print("Excluding outer!")

    for simulation in names:
        filename = "{}/snapshot_{:03d}.hdf5".format(simulation, snapshot)
        if not latex:
            tqdm.write("Attempting to analyse {} {}".format(simulation, snapshot))

        if "lowres" in simulation:
            data = Vertical(filename, res_lr, bbox_x, bbox_y, sf=10, inner=0.5, outer=5)
        else:
            data = Vertical(filename, res_nr, bbox_x, bbox_y, inner=0.5, outer=5)

        name = r"{{\tt {}}} & $0.5<r<5$".format(simulation)
        data.print_data(latex, name)
    
        output_excl[simulation] = data


    
    f, axes = plt.subplots(1, 3)

    for ax, sim in zip(axes, output):
        dat = output[sim]

        ax.hist(dat.dens, bins=100)


    f.savefig("test.png")

""" A large amount of this relies on SurVis as a dependency. It has also been
    stolen from the Jupyter Notebooks in the interstellar-medium-project
    repository also here on GitHub. """

import numpy as np
import matplotlib.pyplot as plt
import survis
from scipy.optimize import curve_fit

# --------------------------- Constants
G = 4.302e-6  # msun kpc km/s
fg = 0.1
f = 0.4
F = 0.5
Pf = 300000  # km/s 100 msun
Qmarg = 1  # marginally stable

#disk
M_disk = 8e9
R_disk = 3

#Nfw
R = 16.1
M_halo = 1.53e12
c = 20.1

#derived
mart_pref = 1.8 * (f/F)**(3/5) * G**(2/5) * Pf**(1/5) * fg**(-2/5)
sigma_0 = M_disk/(2*np.pi*R_disk**2)

# --------------------------------------------------------



def m_nfw(r):
    """ Mass within a radius (NFW profile) -- see Wiki article """

    norm = M_halo/(np.log(1+c) - c/(1+c))

    return (np.log((r + R)/R) - r/(r+R))*norm


def sig_g_custom(r, Q=Qmarg):
    """ Expected Sigma_Gas as a function of radius for a constant Toomre Q
        Parameter Q. """

    this_pref = (np.sqrt(2) * mart_pref)/(np.pi * G * Q  * (1./5.)*(3 + 2./fg))
    v = np.sqrt(G * m_nfw(r)/(r**3))

    return (this_pref * v)**(5/4)


def sig_g_default(r, Q=Qmarg):
    """ Expected Sigma_Gas as a function of radius for a constant toomre Q
        with the default GADGET model """
    v = np.sqrt(G * m_nfw(r)/(r**3))
    sound_speed = np.sqrt((5./3.)*1.38e-23*1e4/1.67e-27)/1e3

    return (sound_speed * np.sqrt(2)* v)/(np.pi * G * Q * (1./5.)*(3 + 2./fg))


def exp_profile(r):
    """ The classic exponential profile """

    return sigma_0 * (np.exp(-r/R_disk))


def Q_sg_r_custom(sg, r):
    """ Used to generate the background colour map, gives Q as a function of
        surface density and radius, for a given galaxy. """

    this_pref = (np.sqrt(2) * mart_pref * 5)/(np.pi * G * (3 + 2./fg) * (sg*1e6)**(4/5))
    v = np.sqrt(G * m_nfw(r)/(r**3))

    return (this_pref * v)

def Q_sg_r_default(sg, r):
    """ Expected Sigma_Gas as a function of radius for a constant toomre Q
        with the default GADGET model """
    v = np.sqrt(G * m_nfw(r)/(r**3))
    sound_speed = np.sqrt((5./3.)*1.38e-23*1e4/1.67e-27)/1e3

    return (sound_speed *np.sqrt(2)* v)/(np.pi * G * (sg*1e6) * (1./5.)*(3 + 2./fg))


def generate_background_matrix(rmin, rmax, sdmin, sdmax, dr, ds, Qfunc=Q_sg_r_custom):
    """ Generates the background matrix using Q_sg_r for the plot's colourmap.
        dr ds give the resolution of the matrix. """

    surface_densities = np.arange(sdmin, sdmax, ds)
    radii = np.arange(rmin, rmax, dr)

    # We need surface_densities as a 2d array and radii as the 'opposite' 2d array

    sd_arr = np.outer(np.flip(surface_densities, 0), np.ones_like(radii))
    rad_arr = np.outer(np.ones_like(surface_densities), radii)

    return Qfunc(sd_arr, rad_arr)


def grab_experimental_data(filename, bins=30):
    """ Grabs the simulation data from file """

    res_elem = 1*(30/bins)
    bbox = [-30, 30]
    res = survis.helper.get_res(res_elem, bbox, bbox)
    data_grid = survis.preprocess.DataGridder(filename, res[0], res[1], bbox[0], bbox[1], bbox[0], bbox[1], False)

    sd_actual = survis.helper.sd_r(data_grid, res_elem, 30, True)

    return np.array([np.array(x) for x in sd_actual])


def get_plottable_exp(filename, bins=30):
    """ SurVis returns things in a pretty dodgy way (sorry about that) so we 
        need to do some creative indexing... """

    raw = grab_experimental_data(filename, bins)

    data = raw[:, 0, 0]
    errors = raw[:, 0, 1]

    return data, errors


def sg_fit(r, Q, x):
    return sig_g(r+x, Q)


def fit_data(filename, bins=30):
    """ Fits the data from (filename) using curve_fit and returns the
        appropriate value for Q and the error. """
    
    data, errors = get_plottable_exp(filename, bins)
    data = data
    errors = errors
    r = np.arange(len(data)*(bins/30))*(30/bins)

    def cut_data(x):
        return x[int(bins/10):int(bins/2)]

    print(cut_data(r))

    popt, perr = curve_fit(sg_fit, cut_data(r), cut_data(data), p0=[0.7, 2], sigma=cut_data(errors), absolute_sigma=True)

    return popt, perr


if __name__ == "__main__":
    # Run in script mode and actually generate the plot
    import matplotlib.gridspec as gridspec
    import matplotlib

    font = {'size':10}
    matplotlib.rc('font', **font)

    np.seterr(invalid='ignore')
    print("Generating the (theory) Toomre Q Figure -- toomre_q.py")

    fig = plt.figure(figsize=(6.3, 3))

    gs = gridspec.GridSpec(1, 3, width_ratios=[10, 10, 1])

    ax_custom = fig.add_subplot(gs[0])
    ax_default = fig.add_subplot(gs[1])
    ax_cb = fig.add_subplot(gs[2])

    rmin = 0
    rmax = 15
    sdmin = 1
    sdmax = 100
    bins = 100
    ds, dr = 0.005, 0.005

    r = np.arange(0, 30, 0.05)
    ax_custom.plot(r, sig_g_custom(r, Q=1)/1e6, color='grey', lw=1)
    ax_default.plot(r, sig_g_default(r, Q=1)/1e6, color='grey', lw=1)

    #ax.plot(r, exp_profile(r)/1e6, label=r"Traditional exponential", ls="dotted")

    # Simulation data
    #gas, err = get_plottable_exp("martizzi_eos_200.hdf5", bins=bins)

    #ax.errorbar(np.arange(len(gas))*(30/bins), gas/1e6, yerr=5*err/1e6, fmt="o", ms=3, label="N-Body Data")
    #popt, pcov = fit_data("martizzi_eos_200.hdf5", bins=bins)
    #ax.plot(r, sg_fit(r, popt[0], popt[1])/1e6, label="Fit to N-Body Data $Q = {:1.3f} \pm {:1.3f}$".format(popt[0], np.sqrt(pcov[0, 0])))

    # Generate the backgorund colormap
    image_data_custom = generate_background_matrix(rmin, rmax, sdmin, sdmax, dr, ds, Q_sg_r_custom)
    im_custom = ax_custom.imshow(image_data_custom,
              extent=[rmin, rmax, sdmin, sdmax],
              vmin=0,
              vmax=2,
              cmap="Spectral",
              aspect="auto")

    image_data_default = generate_background_matrix(rmin, rmax, sdmin, sdmax, dr, ds, Q_sg_r_default)
    im_default = ax_default.imshow(image_data_default,
              extent=[rmin, rmax, sdmin, sdmax],
              vmin=0,
              vmax=2,
              cmap="Spectral",
              aspect="auto")


    for plot in [ax_custom, ax_default]:
        plot.set_ylim(sdmin, sdmax)
        plot.semilogy()
        plot.set_xlim(rmin, rmax)
        plot.set_xlabel("Radius [kpc]")

    ax_default.get_yaxis().set_visible(False)

    ax_custom.set_ylabel("Surface density [$M_\odot$ pc$^{-2}$]")

    plt.colorbar(im_custom, label="Toomre $Q$", cax=ax_cb)

    plt.tight_layout()

    import sys

    if "--showfig" in sys.argv:
        plt.show()
    else:
        plt.savefig("plotgen/toomre_q_theory.pdf", dpi=300)

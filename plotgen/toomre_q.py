""" A large amount of this relies on SurVis as a dependency. It has also been
    stolen from the Jupyter Notebooks in the interstellar-medium-project
    repository also here on GitHub. """

import numpy as np
import matplotlib.pyplot as plt
import survis

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


def sig_g(r, Q=Qmarg):
    """ Expected Sigma_Gas as a function of radius for a constant Toomre Q
        Parameter Q. """

    this_pref = (np.sqrt(2) * mart_pref * fg)/(np.pi * G * Q)
    v = np.sqrt(G * m_nfw(r)/(r**3))

    return (this_pref * v)**(5/4)


def exp_profile(r):
    """ The classic exponential profile """

    return sigma_0 * (np.exp(-r/R_disk))


def Q_sg_r(sg, r):
    """ Used to generate the background colour map, gives Q as a function of
        surface density and radius, for a given galaxy. """

    this_pref = (np.sqrt(2) * mart_pref * fg)/(np.pi * G * (sg*1e6)**(4/5))
    v = np.sqrt(G * m_nfw(r)/(r**3))

    return (this_pref * v)


def generate_background_matrix(rmin, rmax, sdmin, sdmax, dr, ds):
    """ Generates the background matrix using Q_sg_r for the plot's colourmap.
        dr ds give the resolution of the matrix. """

    surface_densities = np.arange(sdmin, sdmax, ds)
    radii = np.arange(rmin, rmax, dr)

    # We need surface_densities as a 2d array and radii as the 'opposite' 2d array

    sd_arr = np.outer(np.flip(surface_densities, 0), np.ones_like(radii))
    rad_arr = np.outer(np.ones_like(surface_densities), radii)

    return Q_sg_r(sd_arr, rad_arr)


if __name__ == "__main__":
    # Run in script mode and actually generate the plot

    fig, ax = plt.subplots()

    rmin = 0
    rmax = 30
    sdmin = 0
    sdmax = 10
    ds, dr = 0.005, 0.01


    r = np.arange(0, 30, 0.05)
    sg_1 = sig_g(r, 1)
    sg_05 = sig_g(r, 0.5)
    sg_5 = sig_g(r, 5)
    ax.plot(r, exp_profile(r)/1e6, label=r"Traditional exponential", ls="dotted")

    im = ax.imshow(generate_background_matrix(rmin, rmax, sdmin, sdmax, dr, ds),
              extent=[rmin, rmax, sdmin, sdmax],
              vmin=0,
              vmax=2,
              cmap="RdYlBu",
              aspect="auto")


    ax.set_ylim(sdmin, sdmax)
    ax.set_xlim(rmin, rmax)
    ax.set_xlabel("Radius [kpc]")
    ax.set_ylabel("Surface density [$M_\odot$ pc$^{-2}$]")
    ax.set_title("Expected $\Sigma_g(r)$ for the Martizzi Model")

    fig.colorbar(im, label="Toomre $Q$")
    
    plt.legend()

    plt.tight_layout()

    plt.show()

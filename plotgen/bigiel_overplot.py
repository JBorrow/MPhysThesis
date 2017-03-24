""" This is the initial test version for the automatic bigiel plot (previously
    I have organised this by hand which is **NOT** ideal """

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

G = 4.3e-3  # pc /msun km^2 /s^2
P = 3000  # km /s
parsecs_in_km = 3.24e-14
seconds_in_year = 31536000
F = 1

prefactor = ((np.pi * G)/(P * F)) * parsecs_in_km * seconds_in_year
logpref = np.log10(prefactor)


def sfr(log_sig_g, fgas):
    return logpref + (2 * log_sig_g) - np.log10(fgas) + 6  # 6 for kpcs


def fgas(sg, tot=50):
    where_too_big = np.greater(sg, tot)
    sg -= sg*where_too_big
    return (sg)/(tot)

if __name__ == "__main__":
    # Run in script mode, actually produce the plot

    fig = plt.figure(figsize=(8,6))

    colors = ['#006388', '#950F2C', '#7E317B', '#164B44', '#C43B8E', '#968E85']

    sgs = np.arange(-1, 5, 0.01)
    
    fg = fgas(10**(sgs))
    sfr_ulirg = sfr(sgs, fg)

    fg100 = fgas(10**(sgs), 100)
    sfr_ulirg100 = sfr(sgs, fg100)

    fg150 = fgas(10**(sgs), 150)
    sfr_ulirg150 = sfr(sgs, fg150)

    fg1000 = fgas(10**(sgs), 1000)
    sfr_ulirg1000 = sfr(sgs, fg1000)

    plt.plot(sgs, sfr(sgs, 0.5), label='$f_\mathrm{gas}$ = 0.5', color=colors[1], ls='dotted', alpha=0.5)
    plt.plot(sgs, sfr(sgs, 0.3), label='$f_\mathrm{gas}$ = 0.3', color=colors[2], ls='dotted', alpha=0.5)
    plt.plot(sgs, sfr(sgs, 0.1), label='$f_\mathrm{gas}$ = 0.1', color=colors[3], ls='dotted', alpha=0.5)

    plt.plot(sgs, sfr_ulirg, color=colors[0], label='$\Sigma_\mathrm{tot}$ = 50 M$_\odot$ pc$^{-2}$')
    plt.plot(sgs, sfr_ulirg100, color=colors[5], label='$\Sigma_\mathrm{tot}$ = 100 M$_\odot$ pc$^{-2}$')
    plt.plot(sgs, sfr_ulirg150, color=colors[4], label='$\Sigma_\mathrm{tot}$ = 150 M$_\odot$ pc$^{-2}$')
    plt.plot(sgs, sfr_ulirg1000, color=colors[4], label='$\Sigma_\mathrm{tot}$ = $10^4$ M$_\odot$ pc$^{-2}$')

    original = Image.open("bigiel_overplot.png")

    plt.imshow(original, extent=[-1, 5, -5, 3], aspect=6/8)

    plt.xlabel(r'log $\Sigma_\mathrm{gas}$ [M$_\odot$ pc$^{-2}$]')
    plt.ylabel(r'log $\dot{\Sigma}_*$ [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]')
    plt.legend(loc=4)

    plt.savefig('bigiel_test.pdf', dpi=600)


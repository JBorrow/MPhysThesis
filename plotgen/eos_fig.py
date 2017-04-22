import numpy as np
import matplotlib.pyplot as plt
import sys

def model(log_nh):
    #springel & hernquist
    if log_nh > -0.89:
        return 0.05*(log_nh**3) - 0.246*(log_nh**2) + 1.749*(log_nh) - 10.6
    else:
        crossing_pressure = 0.05*(-0.89**3) - 0.246*(0.89**2) + 1.749*(-0.89) - 10.6 + 0.89
        return crossing_pressure + log_nh


def EAGLE(nh):
    if nh > 0.1:
        return 4e-12 * (nh**(4./3.))
    else:
        return (1.84e-12) * nh


def my_model(nh):
    if nh > 0.176:
        return 3.43e-12 * (nh**(5./4.))
        #return (6.783e-13/(0.477**(3./2.))) * (nh**(5./4.))
    else:
        return 2.3e-12 * nh


def vertical_line(ax, x=0.1, style='k--', alpha=0.5):
    return ax.plot([x, x], [1e-100, 1e+100], style, alpha=alpha, linewidth=1)


values_of_nh = np.arange(0.01, 100, 0.005)
log_nh = np.log10(values_of_nh)

original = 10**(np.array(list(map(model, log_nh))))
eag = np.array(list(map(EAGLE, values_of_nh)))
my = np.array(list(map(my_model, values_of_nh)))

fig, ax = plt.subplots(figsize=(0.6*6.3, 4))

ax.set_xscale('log')
ax.set_yscale('log')

for line in [0.176, 0.128, 0.1]:
    vertical_line(ax, line)

ax.plot(values_of_nh, original, label="Springel & Hernquist (2003)", linestyle='dotted')
ax.plot(values_of_nh, eag, label="EAGLE, Schaye et. al (2015)", linestyle='-.')
ax.plot(values_of_nh, my, label="This Work")

ax.set_xlabel("$n_H$  [cm$^{-3}$]")
ax.set_ylabel("$P$  [erg cm$^{-3}$]")

ax.set_xlim([0.01, 10])
ax.set_ylim([1e-14, 1e-10])

if "--showfig" in sys.argv:
    fig.show()
    input()
else:
    plt.savefig("plotgen/eos_fig.pdf")

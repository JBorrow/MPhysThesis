""" Creates the SD(r) plot on the server from the pickled data. 
    Unfortunately, because this requires such access to the data,
    it is not possible to make this as part of the executable
    Thesis. """

import matplotlib.pyplot as plt
import numpy as np
import pickle
import survis
from survis.analysis import CommonDataObject


def get_data(filename):
    with open(filename, 'rb') as handle:
        return pickle.load(handle)


def running_average(data, n=2):
    """ Calculates the running average of the data and handles edge cases
        nicely. Returns a numpy array. """

    data = list(data)
    output = []

    # to handle edge cases
    data_edges = data[:n] + data + data[-n:]

    for index in range(len(data)):
        output.append(np.mean(data_edges[index-n:index+n]))

    return np.array(output)


def get_areas(bin_edges):
    output = []

    for i in (np.arange(len(bin_edges)-1)):
        inner = 2*np.pi*bin_edges[i]**2
        outer = 2*np.pi*bin_edges[i+1]**2

        output.append((outer - inner)*1e6)

    return output


def do_plot(ax, data, bin_edges, maxsnap=251, particlemass=80000):
    times = np.arange(maxsnap) * 0.02  # each snap gives 0.02 gyr

    areas = get_areas(bin_edges)

    
    for index, dataset in enumerate(data):
        this_dataset = dataset[:maxsnap] * particlemass
        label = "{} $\leq r <$ {}".format(bin_edges[index], bin_edges[index+1])
        ax.scatter(times, this_dataset/1e9, s=2)
        ax.plot(times, running_average(this_dataset)/1e9, label=label)

    return ax


if __name__ == "__main__":
    # Run in script mode and produce the plot.
    maxsnap = 251
    particlemass = 80000  # msun
    try:
        custom = get_data("custom/processed_variables.pkl")
        default = get_data("default/processed_variables.pkl")
    except FileNotFoundError:
        print("Sorry, but this file is not part of the executable thesis. See the docstring for more information.")
        exit()

    custom_data = survis.analysis.CommonDataExtractor(custom)
    default_data = survis.analysis.CommonDataExtractor(default)

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(6.3, 3))

    do_plot(ax1, custom_data.n_part_r.T, custom_data.bins[0], maxsnap, particlemass)
    do_plot(ax2, default_data.n_part_r.T, default_data.bins[0], maxsnap, particlemass)

    ax1.set_xlim(0, (maxsnap-1)*0.02)
    ax2.set_xlim(0, (maxsnap-1)*0.02)

    ax1.set_ylabel("Mass contained [$10^9$ M$_\odot$]")
    plt.xlabel(r"Time elapsed [Gyr]")

    ax2.legend()
    plt.tight_layout()

    f.savefig("sd_r_fig.pdf")

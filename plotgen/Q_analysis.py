""" Analyses the datai generated by grab_data_many.py on COSMA and produces
    a surface density and Q plot. """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import survis
import pickle

filename = "snapshot_100_data.pkl"

def get_data():
    with open(filename, 'rb') as asset:
        return pickle.load(asset)


def plot_single(ax, data, name, cmap, vmin, vmax, extent, delta=3):
    colormap = cm.get_cmap(cmap)
    c_scale = col.Normalize(vmin=vmin , vmax=vmax)
    colormap.set_bad("white", 1.0)     

    img = ax.imshow(data,
                    vmin=vmin,
                    vmax=vmax,
                    cmap=cmap,
                    extent=extent)

    text = ax.text(extent[0] + delta, extent[2] + delta, name, color="black")

    return ax, img


if __name__ == "__main__":
    """ Run in script mode, and produce plots """
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(6.3, 6.3))
    gs = gridspec.GridSpec(4, 3,
            height_ratios=[10, 10, 1, 9])

    extent = [-30, 30, -30, 30]
    vmin = 0
    vmax = 2
    cmap = "viridis"

    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[0, 2]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[1, 2]),
        fig.add_subplot(gs[2:, 2]),
    ]

    cbar_ax = fig.add_subplot(gs[2, 0:2])

    sim_data = get_data()

    for plot_ax, data in zip(axes, sim_data):
        plot_ax, img = plot_single(plot_ax,
                                   sim_data[data].Q_map,
                                   data,
                                   cmap,
                                   vmin,
                                   vmax,
                                   extent)


    plt.colorbar(img, cax=cbar_ax, orientation='horizontal', label="Toomre $Q$")

    fig.show()
    input()
                                   

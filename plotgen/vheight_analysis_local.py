""" The local analysis of the vertical scale height data generated by
    vertical_scale_height_cluster.py """ 

import numpy as np
import matplotlib.pyplot as plt
import pickle
from vertical_scale_height_cluster import Vertical
from order_of_plots import order

snapshot = 97

def get_data(snapshot=snapshot):
    with open('snapshot_{:03d}_vheight_data.pkl'.format(snapshot), 'rb') as handle:
        return pickle.load(handle)


if __name__ == "__main__":
    # run in script mode

    raw_data = get_data()
    
    fig, ax = plt.subplots()

    for simulation in order:
        this_sim = raw_data[simulation]
        print(this_sim.Z)
        ax.plot(this_sim.r, abs(np.array(this_sim.lZ)), label=simulation)


    plt.legend()

    fig.show()

    input()

        


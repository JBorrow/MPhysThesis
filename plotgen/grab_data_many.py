""" This grabs data from the files on COSMA and makes them into a nice pickled
    format for local development """

import survis
import os
import pickle
from tqdm import tqdm

# ---- Global Parameters
snapshot = 100
bbox_x = [-30, 30]
bbox_y = bbox_x
res_elem = 0.5
res = survis.helper.get_res(res_elem, bbox_x, bbox_y)
# ----------------------


def grab_names():
    """ The simulations live in a directory where they are the only directories
        (but not the only files!) so we must clean os.listdir() """

    output = []
    original = os.listdir()
    
    for filename in original:
        if "." in filename:
            continue
        else:
            output.append(filename)

    return output


if __name__ == "__main__":
    # Run in scripting mode, and actually perform the analysis
    output = {}

    for simulation in tqdm(grab_names()):
        filename = "{}/snapshot_{:03d}.hdf5".format(simulation, snapshot)

        try:
            data = survis.analysis.CommonDataObject(filename, res, bbox_x, bbox_y, res_elem)
        except:
            print("Snapshot {} {} not found, skipping.".format(simulation, snapshot))
            continue

        data.run_analysis()

        output[simulation] = data

    # Now we pickle the data ready for download

    with open('snapshot_{:03d}_data.pkl', 'wb') as pck:
        pickle.dump(result, pck)

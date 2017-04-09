""" Grabs the vertical scale height from the simulations on the cluster and
    saves it to file ready for local analysis. """

import survis
import os
import pickle
from tqdm import tqdm

# ---- Global Parameters
snapshot = 97
bbox_x = [-40, 40]
bbox_y = bbox_x
res_elem_lr = 0.5
res_elem_nr = 0.2
res_lr = survis.helper.get_res(res_elem_lr, bbox_x, bbox_y)
res_nr = survis.helper.get_res(res_elem_nr, bbox_x, bbox_y)
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



class Vertical(object):
    def __init__(self, filename, res, bbox_x, bbox_y, max_rad=30, n_bins=30):
        data = survis.preprocess.DataGridder(filename,
                                             res[0],
                                             res[1],
                                             bbox_x[0],
                                             bbox_x[1],
                                             bboy_y[0],
                                             bboy_y[1],
                                             False)

        self.filename = filename
        self.res = res
        self.bbox_x = bbox_x
        self.bbox_y = bbox_y
        self.max_rad = max_rad
        self.n_bins = n_bins

        crds = data['PartType0']['Coordinates']
        z_data = data['PartType0']['Coordinates'][:, 2]

        radii = np.sqrt(np.sum(np.square(crds), 1))
        
        self.r, self.Z, self.dZ = self.scale_height_r(z_data, radii, n_bins, max_rad)

        return


    def bin_cent(self, bins):
        return 0.5*(bins[1:] + bins[:1])


    def get_scale_height(self, z_data, min=-0.4, max=0.4, nbins=50):

        def to_fit(z, norm, Z):
            return norm*(1/(np.cosh(z/Z)**2))

        # First we must bin the data
        bins = np.arange(min, max, (max-min)/nbins)
        n, bins = np.arange(z_data, bins)
        bincenters = bin_cent(bins)

        popt, pcov = curve_fit(to_fit, bincenters, n)

        return popt[1], np.sqrt(pcov[1, 1])


    def scale_height_r(self, z_data, radii, n_bins, max_rad):
        """ Calculates the vertical scale height as a function of radius in the
        disk by binning the data twice """

        radii_bins = np.arange(0, max_rad, max_rad/n_bins)
        output_Z = []
        output_dZ = []

        for index in range(len(radii_bins) - 1):  # We want to skip the first and the last
            bottom = radii_bins[index+1]
            top = radii_bins[index+2]

            mask = np.logical_or((radii < bottom), (radii > top))
            current_data = np.ma.masked_array(z_data, mask).compressed()

            Z, dZ = self.get_scale_height(current_data)

            output_Z.append(Z)
            output_dZ.append(dZ)

        return self.bin_centers(radii_bins), output_Z, output_dZ


if __name__ == "__main__":
    # Run in scripting mode, and actually perform the analysis
    output = {}

    for simulation in tqdm(grab_names()):
        filename = "{}/snapshot_{:03d}.hdf5".format(simulation, snapshot)
        tqdm.write("Attempting to analyse {} {}".format(simulation, snapshot))

        if "lowres" in simulation:
            data = Vertical(filename, res_lr, bbox_x, bbox_y)
        else:
            data = Vertical(filename, res_hr, bbox_x, bbox_y)
        tqdm.write("Analysis of {} {} successful.".format(simulation, snapshot))
    
        output[simulation] = data

    # Now we pickle the data ready for download

    with open('snapshot_{:03d}_vheight_data.pkl'.format(snapshot), 'wb') as pck:
        pickle.dump(output, pck)


import numpy as np


def get_coords(in_file):
    '''
    Return coordinates from file.
    '''
    if in_file == 'random.dat':
        # Generate random data.
        N = 5000
        ra = np.random.uniform(low=0., high=1000., size=(N))
        dec = np.random.uniform(low=0., high=1000., size=(N))
        mag = np.random.uniform(low=10., high=24., size=(N))
    else:
        # Each sub-list in 'in_file' is a row of the file.
        file_data = np.loadtxt(in_file)

        # Extract coordinates and zip them into two lists.
        ra, dec, mag = zip(*file_data)[1], zip(*file_data)[2], \
            zip(*file_data)[3]

    return ra, dec, mag


import numpy as np


def get_coords(in_file, mag_range):
    '''
    Return coordinates from file.
    '''
    if in_file == 'random.dat':
        # Generate random data.
        N = 5000
        x = np.random.uniform(low=0., high=1000., size=(N))
        y = np.random.uniform(low=0., high=1000., size=(N))
        mag = np.random.uniform(low=10., high=24., size=(N))
    else:
        # Each sub-list in 'in_file' is a row of the file.
        file_data = np.loadtxt(in_file)

        # Extract coordinates and zip them into lists.
        x, y, mag = zip(*file_data)[1], zip(*file_data)[2], \
            zip(*file_data)[3]

    # Filter stars outside of magnitude range
    x_mr, y_mr = [], []
    for i, m in enumerate(mag):
        if mag_range[0] <= m < mag_range[1]:
            x_mr.append(x[i])
            y_mr.append(y[i])

    return x_mr, y_mr, x, y, mag

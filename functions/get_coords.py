
import numpy as np


def get_coords(in_file, in_file_cols, mag_range, coords_flag):
    '''
    Return coordinates from file.
    '''
    if in_file == 'random.dat':
        # Generate random data.
        N = 10000
        x = np.random.uniform(low=0., high=1000., size=(N))
        y = np.random.uniform(low=0., high=1000., size=(N))
        mag = np.random.uniform(low=10., high=24., size=(N))
    else:
        # Each sub-list in 'in_file' is a row of the file.
        file_data = np.loadtxt(in_file)
        i, j, k = in_file_cols

        # Extract coordinates and zip them into lists.
        x, y, mag = zip(*file_data)[i], zip(*file_data)[j], \
            zip(*file_data)[k]

    if coords_flag == 'deg':
        ra_cent = (max(x) + min(x)) / 2.
        dec_cent = (max(y) + min(y)) / 2.
    else:
        ra_cent, dec_cent = 0., 0.

    # Move coordinates to origin in center of frame. Apply cos(Dec) correction
    # to RA coordinates.
    x = (np.array(x) - ra_cent) * np.cos(np.deg2rad(dec_cent))
    y = np.array(y) - dec_cent

    # Filter stars outside of magnitude range
    x_mr, y_mr, mag_mr = [], [], []
    for i, m in enumerate(mag):
        if mag_range[0] <= m < mag_range[1]:
            x_mr.append(x[i])
            y_mr.append(y[i])
            mag_mr.append(m)

    return x_mr, y_mr, mag_mr, x, y, mag, ra_cent, dec_cent

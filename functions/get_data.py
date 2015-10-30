
import numpy as np


def get_cents_rad_file(coords_flag, vor_flag, cr_file_cols, ra_cent, dec_cent):
    '''
    Load center coordinates and radii values from file.
    '''
    # Load file with centers and radii.
    com_file = vor_flag
    file_data = np.loadtxt(com_file)
    i, j, k, q = cr_file_cols
    # Extract coordinates and zip them into lists.
    cc_x, cc_y, cr_x, cr_y = zip(*file_data)[i], zip(*file_data)[j], \
        zip(*file_data)[k], zip(*file_data)[q]

    if coords_flag != 'deg':
        ra_cent, dec_cent = 0., 0.
        arcmin_2_deg = 1.
    else:
        arcmin_2_deg = 60. * 2.

    # Move center coordinates to new origin.
    cc_x = (np.array(cc_x) - ra_cent) * np.cos(np.deg2rad(dec_cent))
    cc_y = (np.array(cc_y) - dec_cent)

    # Convert RA, DEC values total box longitudes in arcmin to
    # radii in degrees.
    cr_x, cr_y = np.array(cr_x) / arcmin_2_deg, np.array(cr_y) / arcmin_2_deg

    # Select larger radius value.
    cc_r = []
    for rx, ry in zip(*[cr_x, cr_y]):
        cc_r.append(max(rx, ry))
    # Zip data.
    cent_rad = zip(*[cc_x, cc_y, cc_r])

    return cent_rad


def get_data(in_file, in_file_cols, mag_range, coords_flag, vor_flag):
    '''
    Return coordinates from file.
    '''
    if in_file == 'random.dat':
        # Generate random data.
        N = 10000
        x = np.random.uniform(low=0., high=1000., size=(N))
        y = np.random.uniform(low=0., high=1000., size=(N))
        mag = np.random.uniform(low=10., high=22., size=(N))
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


import numpy as np
from math import hypot


def in_radius(c_x, c_y, r, x, y):
    return hypot(c_x - x, c_y - y) <= r


def filt_density(fr_area, x_mr, y_mr, mag_mr, cent_rad, dens_frac):
    '''
    Apply density filter.
    '''

    # Frame's density.
    fr_dens_area = len(x_mr) / fr_area

    # Obtain integrated magnitude for each defined circle.
    dens_accp_groups, dens_rej_groups = [], []
    for c_x, c_y, r in cent_rad:
        # Count stars within this circle, using stars that passed the
        # magnitude filter.
        clust_mags, N = [], 0
        for i, (x, y) in enumerate(zip(*[x_mr, y_mr])):
            if in_radius(c_x, c_y, r, x, y):
                # This is used later on in the I/A filter function.
                clust_mags.append(mag_mr[i])
                N += 1

        # Obtain density for this cluster, normalized to 1. for the frame's
        # density.
        clust_dens = (float(N) / (np.pi * (r ** 2))) / fr_dens_area

        # If the overdensity has a density larger than a given fraction of
        # the frame's density, keep (j=0). Else, discard (j=1).
        if clust_dens > dens_frac:
            dens_accp_groups.append([c_x, c_y, r, clust_dens, clust_mags])
        else:
            dens_rej_groups.append([c_x, c_y, r, clust_dens, clust_mags])

    return dens_accp_groups, dens_rej_groups

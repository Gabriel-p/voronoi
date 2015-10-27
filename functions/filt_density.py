
import numpy as np
from math import hypot


def in_radius(c_x, c_y, r, x, y):
    return hypot(c_x - x, c_y - y) <= r


def filt_density(fr_area, x_mr, y_mr, cent_rad, dens_frac):
    '''
    Apply density filter.
    '''

    # Frame's density.
    fr_dens_area = len(x_mr) / fr_area

    # Obtain integrated magnitude for each defined circle.
    dens_accp_groups, dens_rej_groups, dens_all = [], [], [[], [], []]
    for c_x, c_y, r in cent_rad:
        # Count stars within this circle, using stars that passed the
        # magnitude filter.
        N = 0
        for x, y in zip(*[x_mr, y_mr]):
            if in_radius(c_x, c_y, r, x, y):
                N += 1

        # Obtain density for this cluster, normalized to 1. for the frame's
        # density.
        clust_dens = (float(N) / (np.pi * (r ** 2))) / fr_dens_area

        # If the overdensity has a density larger than a given fraction of
        # the frame's density, keep (j=0). Else, discard (j=1).
        if clust_dens > dens_frac:
            dens_accp_groups.append([c_x, c_y, r])
            dens_all[0].append(clust_dens)
        else:
            dens_rej_groups.append([c_x, c_y, r])
            dens_all[1].append(clust_dens)
        # Store all values for plotting.
        dens_all[2].append(clust_dens)

    return dens_accp_groups, dens_rej_groups, dens_all

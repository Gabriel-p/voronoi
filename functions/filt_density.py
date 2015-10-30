
import numpy as np
from math import hypot
from save_to_file import save_to_log


def in_radius(c_x, c_y, r, x, y):
    return hypot(c_x - x, c_y - y) <= r


def filt_density(f_name, fr_area, x_mr, y_mr, mag_mr, cent_rad, vor_flag,
                 area_frac_range, dens_frac, mean_filt_area):
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

    # If the center coords and radii were read from file, obtain the min and
    # max fraction of mean Voronoi area for these clusters.
    if vor_flag != 'voronoi':
        all_fracs = []
        # Obtain, for each accepted cluster, its area/N divided by the mean
        # Voronoi area value.
        for g in dens_accp_groups:
            cl_dens = (1. / fr_dens_area) * (1. / g[3])
            all_fracs.append(cl_dens / mean_filt_area)
        area_frac_range = [min(all_fracs), max(all_fracs)]
        save_to_log(f_name, 'Obtained area range: [{:.2f}, {:.2f}]'.format(
            *area_frac_range), 'a')

    return dens_accp_groups, dens_rej_groups, area_frac_range

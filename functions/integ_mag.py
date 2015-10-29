
import numpy as np


def calc_integ_mag(mags):
    '''
    Calculate integrated magnitude up to a certain maximum magnitude value.
    '''

    # Sort magnitude list.
    sort_lis = sorted(mags)

    int_mag_val = 0.

    energ_sum = 0.
    for mag_i in sort_lis:
        energ_sum = energ_sum + 10 ** (mag_i / -2.5)

    int_mag_val = -2.5 * np.log10(max(1e-20, energ_sum))

    return int_mag_val


def filt_integ_mag(x_mr, y_mr, mag_mr, dens_accp_groups, dens_rej_groups,
                   intens_frac):
    '''
    Apply integrated magnitude filter.
    '''
    # Obtain frame's integrated magnitude per area unit.
    frame_int_mag = calc_integ_mag(mag_mr)
    frame_area = (max(x_mr) - min(x_mr)) * (max(y_mr) - min(y_mr))
    frame_intens_area = 1.

    # Obtain integrated magnitude for each defined circle.
    intens_dens_all = [[[[], [], [], []], [[], [], [], []]],
                       [[[], [], [], []], [[], [], [], []]]]
    for i, ls in enumerate([dens_accp_groups, dens_rej_groups]):
        for (c_x, c_y, r, clust_dens, clust_mags) in ls:
            # Integrated magnitude.
            clust_int_mag = calc_integ_mag(clust_mags)
            # Area
            clust_area = np.pi * (r ** 2)
            # Intensity
            clust_intens = 10 ** (0.4 * (frame_int_mag - clust_int_mag))
            # I/A
            clust_intens_area = (clust_intens * frame_area) / clust_area

            # If the overdensity has an intensity per unit area (I/A) larger
            # than a given fraction of the frame's I/A keep, else, reject.
            if clust_intens_area > intens_frac * frame_intens_area:
                intens_dens_all[i][0][0].append([c_x, c_y, r])
                intens_dens_all[i][0][1].append(len(clust_mags))
                intens_dens_all[i][0][2].append(clust_dens)
                intens_dens_all[i][0][3].append(clust_intens_area)
            else:
                intens_dens_all[i][1][0].append([c_x, c_y, r])
                intens_dens_all[i][1][1].append(len(clust_mags))
                intens_dens_all[i][1][2].append(clust_dens)
                intens_dens_all[i][1][3].append(clust_intens_area)

    # Pack into separated lists.
    intens_acc_dens_acc, intens_acc_dens_rej, intens_rej_dens_acc,\
        intens_rej_dens_rej = intens_dens_all[0][0], intens_dens_all[1][0],\
        intens_dens_all[0][1], intens_dens_all[1][1]

    return intens_acc_dens_acc, intens_acc_dens_rej, intens_rej_dens_acc,\
        intens_rej_dens_rej


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


def filt_integ_mag(x_mr, y_mr, pts_area_thres, mag_area_thres, cent_rad,
                   intens_frac):
    '''
    Apply integrated magnitude filter.
    '''
    # Unpack
    x_thr, y_thr = zip(*pts_area_thres)

    # Store *all* intensities per area unit.
    intens_area_all = [[[], [], []], [[], [], []], []]

    # Obtain frame's integrated magnitude per area unit.
    frame_int_mag = calc_integ_mag(mag_area_thres)
    frame_area = (max(x_thr) - min(x_thr)) * (max(y_thr) - min(y_thr))
    frame_intens_area = 1.
    # print "Frame's integrated magnitude: {}".format(frame_int_mag)

    # Obtain integrated magnitude for each defined circle.
    intes_accp_groups, intens_rej_groups, clust_intens_area = [], [], []
    for c_x, c_y, r in cent_rad:
        # Count stars within this circle, using stars that passed the
        # magnitude filter.
        N = 0
        for x, y in zip(*[x_mr, y_mr]):
            d = np.sqrt((c_x - x) ** 2 + (c_y - y) ** 2)
            if d <= r:
                N += 1

        # Group stars within this circle, using only stars that passed the area
        # threshold.
        clust_mags = []
        for i, (x, y) in enumerate(pts_area_thres):
            d = np.sqrt((c_x - x) ** 2 + (c_y - y) ** 2)
            if d <= r:
                clust_mags.append(mag_area_thres[i])
        clust_int_mag = calc_integ_mag(clust_mags)
        clust_area = np.pi * (r ** 2)
        clust_intens = 10 ** (0.4 * (frame_int_mag - clust_int_mag))
        clust_intens_area = (clust_intens * frame_area) / clust_area

        # If the overdensity has an intensity per unit area (I/A) larger than a
        # given fraction of the frame's I/A, keep (j=0). Else, discard (j=1).
        if clust_intens_area > intens_frac * frame_intens_area:
            intes_accp_groups.append([c_x, c_y, r])
            j = 0
        else:
            intens_rej_groups.append([c_x, c_y, r])
            j = 1
        # Store intensity/area for plotting.
        intens_area_all[j][0].append(clust_intens_area)
        intens_area_all[j][1].append(r)
        intens_area_all[j][2].append(N)
        # Store all values for plotting.
        intens_area_all[2].append(clust_intens_area)

    return intes_accp_groups, intens_rej_groups, intens_area_all

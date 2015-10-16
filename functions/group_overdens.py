
import numpy as np
from itertools import combinations
import operator


def merge_overdens(cent_rad):
    '''
    Take the groups found and compare each one with all the others.

    If the distance between the centers of two groups is less than the maximum
    radius value of the two, the groups are considered the same one.

    Keep the group with the larger radius.
    '''

    # sort list of center and radii putting those with smaller radius first.
    sort_cent_rad = sorted(cent_rad, key=operator.itemgetter(2))

    keep_idx, disc_idx, keep_flag = [], [], True
    # Store lists and their indexes.
    l2 = [(i, sort_cent_rad[i]) for i in range(len(sort_cent_rad))]
    # Loop through every combination of groups/overdensities/elements.
    for (i, a), (j, b) in combinations(l2, 2):
        # Unpack values.
        c_x, c_y, r = a  # 'i' group values
        n_cx, n_cy, n_r = b  # 'j' group values

        # Distance between the two groups.
        d = np.sqrt((n_cx - c_x) ** 2 + (n_cy - c_y) ** 2)
        # If the two centers are closer than the maximum radius.
        if d <= max(r, n_r):
            # If the new 'j' group has a larger radius.
            if n_r > r:
                # Don't keep the current overdensity.
                keep_flag = False

        # Reached the end of group's 'i' comparison with the remaining
        # elements.
        if j + 1 == len(sort_cent_rad):
            # If no group was found close to this one and with a larger radius,
            # keep the current group.
            if keep_flag:
                keep_idx.append(i)
            else:
                # Discard this group.
                disc_idx.append(i)
            # Reset flag.
            keep_flag = True

    # Store index of last item (overdensity with the largest radius)
    keep_idx.append(len(sort_cent_rad) - 1)

    # Store groups that were kept.
    new_cent_rad = [sort_cent_rad[k] for k in keep_idx]
    # Store groups that were merged and discarded, for plotting.
    old_cent_rad = [sort_cent_rad[i] for i in disc_idx]

    return old_cent_rad, new_cent_rad

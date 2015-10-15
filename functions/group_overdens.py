
import numpy as np
from itertools import combinations


# def check_cent_rad(cent_rad, cent_rad_all):
#     '''
#     Take new groups found with a new minimum neighbors value, compare with
#     existing groups (obtained with larger m_n values) and decide which one to
#     keep.

#     Criteria: if the distance between the centers is less than the maximum
#     radius value, then the two groups are considered the same one.
#     If the old group has a larger radius, then that one is kept. If the new
#     group has a larger radius, then the new one is kept.
#     '''

#     new_cent_rads = []
#     # For each new overdensity identified.
#     for n_cx, n_cy, n_r in cent_rad:
#         new_overdens, overwrite_overdens = False, False

#         # For each already stored overdensity.
#         for idx, [c_x, c_y, r] in enumerate(cent_rad_all):

#             # Only process if the groups are not the same.
#             if [n_cx, n_cy, n_r] != [c_x, c_y, r]:
#                 d = np.sqrt((n_cx - c_x) ** 2 + (n_cy - c_y) ** 2)
#                 # See if the center of one overdensity is within the area of
#                 # the other.
#                 if d <= max(r, n_r):
#                     # If the old radius is larger than the new one.
#                     if r >= n_r:
#                         # Keep the old overdensity.
#                         print 'Keep (x={:.2f}, y={:.2f}, r={:.2f})'.format(
#                             c_x, c_y, r)
#                         print ' versus (x={:.2f}, y={:.2f}, r={:.2f})'.format(
#                             n_cx, n_cy, n_r)
#                     else:
#                         # Keep the new one and discard the old.
#                         ow_idx = idx
#                         overwrite_overdens = True
#                         print 'Replace (x={:.2f}, y={:.2f}, r={:.2f})'.format(
#                             c_x, c_y, r)
#                         print ' with (x={:.2f}, y={:.2f}, r={:.2f})'.format(
#                             n_cx, n_cy, n_r)

#         # If the new overdensity was not found within any known one, add it.
#         if new_overdens:
#             if overwrite_overdens:
#                 cent_rad_all[ow_idx] = [n_cx, n_cy, n_r]
#             else:
#                 # No match with a known overdensity. Add to list.
#                 # print 'Add (cx={:.2f}, cy={:.2f}, r={:.2f})'.format(
#                 #     n_cx, n_cy, n_r)
#                 new_cent_rads.append([n_cx, n_cy, n_r])

#     return new_cent_rads, cent_rad_all


def merge_overdens(cent_rad):
    '''
    Take the groups found with a fixed minimum neighbors value and compare
    each one with all the others.

    If the distance between the centers of two groups is less than the maximum
    radius value of the two, the groups are considered the same one.

    Keep the group with the larger radius.
    '''
    # Generate a copy of the list.
    new_cent_rad = [_ for _ in cent_rad]

    keep_grouping = True
    while keep_grouping:
        keep_grouping = False

        keep_idx, idx = [], 'none'
        # Store lists and their indexes.
        l2 = [(i, new_cent_rad[i]) for i in range(len(new_cent_rad))]
        for (i, a), (j, b) in combinations(l2, 2):
            # Unpack values.
            c_x, c_y, r = a
            n_cx, n_cy, n_r = b

            # Distance between the two groups.
            d = np.sqrt((n_cx - c_x) ** 2 + (n_cy - c_y) ** 2)
            # If the two centers are closer than the maximum radius.
            if d <= max(r, n_r):
                # If the current radius is larger than the previous one.
                if n_r > r:
                    # Keep the index of the current overdensity.
                    idx = j

            # Reached the end of group's i comparison with the remaining
            # elements.
            if j + 1 == len(new_cent_rad):
                # If some group was found within the minimum distance and with
                # a larger radius, store its index.
                idx_s = idx if idx != 'none' else i
                # If this index was not stored previously, store it.
                if idx_s not in keep_idx:
                    keep_idx.append(idx_s)
                # Reset index to be stored.
                idx = 'none'

        # Store index of last item, if it was not stored previously.
        if len(new_cent_rad) - 1 not in keep_idx:
            keep_idx.append(len(new_cent_rad) - 1)

        # If at least 1 group was merged, iterate merging algorithm again.
        if len(keep_idx) < len(new_cent_rad):
            keep_grouping = True
            new_cent_rad = [new_cent_rad[k] for k in keep_idx]

    # Store groups that were merged and discarded, for plotting.
    old_cent_rad = [i for i, j in zip(cent_rad, new_cent_rad) if i != j]

    return old_cent_rad, new_cent_rad

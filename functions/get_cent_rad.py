
import numpy as np
from smallestenclosingcircle import make_circle
from get_data import get_cents_rad_file
from save_to_file import save_to_log
from percent_done import print_perc


def area_filter(acp_pts, acp_mags, pts_area, pts_vert, mean_filt_area,
                area_frac_range):
    '''
    Filter *out* those points whose area is *outside a range given by a
    fraction of the mean Voronoi cells area value passed.

    http://stackoverflow.com/a/32058576/1391441
    '''

    pts_thres, mag_thres, vert_thres = [], [], []
    for i, p in enumerate(acp_pts):
        # Keep point if its area is within the mgiven threshold.
        if mean_filt_area * area_frac_range[0] < pts_area[i] < \
                mean_filt_area * area_frac_range[1]:
            pts_thres.append(p)
            mag_thres.append(acp_mags[i])
            vert_thres.append(pts_vert[i])

    return pts_thres, mag_thres, vert_thres


def consolidate(sets):
    # http://rosettacode.org/wiki/Set_consolidation#Python:_Iterative
    setlist = [s for s in sets if s]

    tot_sols = len(setlist)
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    for i, s1 in enumerate(setlist):
        if s1:
            for s2 in setlist[i+1:]:
                intersection = s1.intersection(s2)
                if intersection:
                    s2.update(s1)
                    s1.clear()
                    s1 = s2

        # Print percentage done.
        milestones = print_perc(i, tot_sols, milestones)

    return [s for s in setlist if s]


def wrapper(seqs):
    consolidated = consolidate(map(set, seqs))
    groupmap = {x: i for i, seq in enumerate(consolidated) for x in seq}

    output = {}
    for seq in seqs:
        target = output.setdefault(groupmap[seq[0]], [])
        target.append(seq)

    return list(output.values())


def shared_vertex(vert_thres):
    '''
    Check which of the points that passed the area filter share at least
    one vertex. If they do, they are considered neighbors and stored in the
    same list.
    '''
    all_groups = wrapper(vert_thres)
    return all_groups


def element_count(p1, p2):
    '''
    Count points in nested lists.
    '''
    count = 0
    # For each group defined.
    for g in p1:
        # For each point in each group.
        for _ in g:
            # For each point above threshold.
            count += len(p2)

    return count


def group_stars(pts_thres, mag_thres, vert_thres, all_groups):
    '''
    For each defined group, find the points that correspond to it.
    '''

    tot_sols = element_count(all_groups, pts_thres)
    milestones, c = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 0

    neighbors = [[] for _ in range(len(all_groups))]
    neig_mags = [[] for _ in range(len(all_groups))]
    # For each group defined.
    for i, g in enumerate(all_groups):
        # For each point in each group.
        for vert_pt in g:
            # For each point above threshold.
            for j, p in enumerate(pts_thres):
                c += 1
                verts = vert_thres[j]
                if verts == vert_pt:
                    neighbors[i].append(p)
                    neig_mags[i].append(mag_thres[j])

                # Print percentage done.
                milestones = print_perc(c, tot_sols, milestones)

    return neighbors, neig_mags


def neighbors_filter(neighbors, neig_mags, m_m_n):
    '''
    Filter the points that passed the area filter, discarding those which
    have fewer neighbors than the 'min_neighbors' value.
    '''
    pts_neighbors, mags_neighbors = [], []
    for i, g in enumerate(neighbors):
        if m_m_n[0] <= len(g) <= m_m_n[1]:
            pts_neighbors.append(zip(*g))
            for m in neig_mags[i]:
                mags_neighbors.append(m)

    return pts_neighbors, mags_neighbors


def get_cent_rad(f_name, coords_flag, cr_file_cols, m_m_n, vor_flag, ra_cent,
                 dec_cent, acp_pts, acp_mags, pts_area, pts_vert,
                 mean_filt_area, area_frac_range):
    '''
    Assign a center and a radius for each overdensity/group identified.

    Use an automatic algorithm based on a Voronoi diagram and grouping neighbor
    stars, or a file with center coordinates and radii already stored.
    '''

    if vor_flag == 'voronoi':

        # Apply maximum area filter.
        pts_area_thres, mag_area_thres, vert_area_thres = area_filter(
            acp_pts, acp_mags, pts_area, pts_vert, mean_filt_area,
            area_frac_range)
        save_to_log(f_name, "\nStars filtered by area range: {}".format(
            len(pts_area_thres)), 'a')

        save_to_log(f_name, 'Detect shared vertex', 'a')
        all_groups = shared_vertex(vert_area_thres)

        # This is a costly process.
        save_to_log(f_name, 'Assign points to groups', 'a')
        neighbors, neig_mags = group_stars(pts_area_thres, mag_area_thres,
                                           vert_area_thres, all_groups)
        save_to_log(f_name, 'Groups found: {}'.format(len(neighbors)), 'a')

        # Keep only those groups with a higher number of members than
        # min_neighbors.
        pts_neighbors, mags_neighbors = neighbors_filter(neighbors, neig_mags,
                                                         m_m_n)
        save_to_log(f_name,
                    '\nGroups filtered by number of members: {}'.format(
                        len(pts_neighbors)), 'a')

        # Check if at least one group was defined with the minimum
        # required number of neighbors.
        if pts_neighbors:
            # Obtain center and radius for each overdensity identified.
            cent_rad = []
            for g in pts_neighbors:

                pts = np.array(zip(*g))
                mean_pts = np.mean(pts, 0)

                # Translate towards origin
                pts -= mean_pts
                result = make_circle(pts)

                # Move back to correct position.
                c_r = (result[0] + mean_pts[0], result[1] + mean_pts[1],
                       result[2])

                cent_rad.append([c_r[0], c_r[1], c_r[2]])
        else:
            cent_rad = []
            save_to_log(f_name, "No groups found with members in range {}".
                        format(m_m_n), 'a')

    else:
        # Get data from file.
        cent_rad = get_cents_rad_file(coords_flag, vor_flag, cr_file_cols,
                                      ra_cent, dec_cent)
        pts_area_thres, mag_area_thres = [[0., 0.], [0., 0.]], [0.]
        pts_neighbors, mags_neighbors = [[[0.], [0.]]], [0]

    return cent_rad, pts_area_thres, mag_area_thres, pts_neighbors,\
        mags_neighbors

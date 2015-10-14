

from scipy.spatial import Voronoi
import time
import numpy as np
from functions.get_params_in import get_params_in
from functions.get_coords import get_coords
from functions.get_cent_rad import get_cent_rad
from functions.save_to_file import save_cent_rad, save_to_log
from functions.make_plots import area_hist, vor_plot


def print_perc(i, tot_sols, milestones):
    '''
    '''
    percentage_complete = (100.0 * (i + 1) / tot_sols)
    while len(milestones) > 0 and \
            percentage_complete >= milestones[0]:
        print " {:>3}% done".format(milestones[0])
        # Remove that milestone from the list.
        milestones = milestones[1:]

    return milestones


def area_of_polygon(points):
    """
    Calculates the area of an arbitrary polygon given its vertices

    Source: http://stackoverflow.com/a/4682656/1391441
    """
    x, y = zip(*points)
    area = 0.0
    for i in xrange(-1, len(x)-1):
        area += x[i] * (y[i+1] - y[i-1])
    return abs(area) / 2.0


def get_vor_data(points, vor):
    '''
    Obtain from Voronoi diagram: acp_pts (points whose diagram/cell/polygon
    contains no negative indexes), rej_pts (points with negative index vertex),
    pts_area (area of each cell associated with each accepted point),
    pts_vert (indexes of vertices associated with each point).
    '''

    tot_sols = len(vor.regions)
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    acp_pts, rej_pts, pts_area, pts_vert = [], [], [], []
    for i, reg in enumerate(vor.regions):
        # If region is not empty.
        if reg:

            # Get index of the point that corresponds to this region.
            p_idx = np.where(vor.point_region == i)[0][0]

            # Check that the vertices index for the region are all positive.
            if all([_ >= 0 for _ in reg]):
                # Store point coordinates.
                acp_pts.append(points[p_idx])

                # Store coordinates for all vertices.
                pts = []
                for p in reg:
                    vert = vor.vertices[p]
                    pts.append(vert)

                # Get cell/region/polygon area.
                cell_area = area_of_polygon(pts)
                pts_area.append(cell_area)

                # Store region's vertex index.
                pts_vert.append(reg)
            else:
                rej_pts.append(points[p_idx])

        # Print percentage done.
        milestones = print_perc(i, tot_sols, milestones)

    return acp_pts, rej_pts, pts_area, pts_vert


def area_filter(acp_pts, pts_area, pts_vert, avr_area, avr_area_frac):
    '''
    Filter those points whose area is below a certain fraction of the average
    area.

    http://stackoverflow.com/a/32058576/1391441
    '''

    pts_thres, vert_thres = [], []
    for i, p in enumerate(acp_pts):
        if pts_area[i] < avr_area_frac * avr_area:
            pts_thres.append(p)
            vert_thres.append(pts_vert[i])

    return pts_thres, vert_thres


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


def group_stars(pts_thres, vert_thres, all_groups):
    '''
    For each defined group, find the points that correspond to it.
    '''

    tot_sols = len(all_groups)
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    neighbors = [[] for _ in range(len(all_groups))]
    # For each group defined.
    for i, g in enumerate(all_groups):
        # For each point in each group.
        for vert_pt in g:
            # For each point above threshold.
            for j, p in enumerate(pts_thres):
                verts = vert_thres[j]
                if verts == vert_pt:
                    neighbors[i].append(p)

        # Print percentage done.
        milestones = print_perc(i, tot_sols, milestones)

    return neighbors


def neighbors_filter(neighbors, min_neighbors):
    '''
    Filter the points that passed the area filter, discarding those which
    have fewer neighbors than the 'min_neighbors' value.
    '''

    # Keep only those groups with a higher number of members than
    # min_neighbors.
    pts_neighbors = []
    for g in neighbors:
        if len(g) > min_neighbors:
            pts_neighbors.append(zip(*g))

    return pts_neighbors


def check_cent_rad(cent_rad, cent_rad_all):
    '''
    Take new groups found with a new minimum neighbors value, compare with
    existing groups (obtained with larger m_n values) and decide which one to
    keep.

    Criteria: if the distance between the centers is less than the maximum
    radius value, then the two groups are considered the same one.
    If the old group has a larger radius, then that one is kept. If the new
    group has a larger radius, then the new one is kept.
    '''

    new_cent_rads = []
    # for each new overdensity identified.
    for n_cx, n_cy, n_r in cent_rad:
        new_overdens, overwrite_overdens = True, False
        # for each already stored overdensity.
        for idx, [c_x, c_y, r] in enumerate(cent_rad_all):
            d = np.sqrt((n_cx - c_x) ** 2 + (n_cy - c_y) ** 2)
            # See if the new overdensity is located inside this already
            # found overdensity.
            # The center of one overdensity is within the area of the other.
            if d <= max(r, n_r):
                # If the old radius is larger than the new one.
                if r >= n_r:
                    # Keep the old overdensity.
                    new_overdens = False
                    print 'Keep (cx={:.2f}, cy={:.2f}, r={:.2f})'.format(
                        c_x, c_y, r)
                else:
                    # Keep the new one and discard the old.
                    ow_idx = idx
                    overwrite_overdens = True
                    print 'Replace (cx={:.2f}, cy={:.2f}, r={:.2f})'.format(
                        c_x, c_y, r)
                    print 'with (cx={:.2f}, cy={:.2f}, r={:.2f})'.format(
                        n_cx, n_cy, n_r)
            else:
                # No match with this known overdensity. Add to list.
                print 'Add (cx={:.2f}, cy={:.2f}, r={:.2f})'.format(
                    n_cx, n_cy, n_r)
                pass

        # If the new overdensity was not found within any known one, add it.
        if new_overdens:
            if overwrite_overdens:
                cent_rad_all[ow_idx] = [n_cx, n_cy, n_r]
            else:
                new_cent_rads.append([n_cx, n_cy, n_r])

    return new_cent_rads, cent_rad_all


def main():
    # Read parameters from input file.
    in_file, mag_range, avr_area_frac, min_neighbors = get_params_in()

    # Each sub-list in 'in_file' is a row of the file.
    f_name = in_file[:-4]
    # Create log file.
    text = 'Processing file: {}\n'.format(f_name)
    print text
    save_to_log(f_name, text, 0)
    text = ''

    # Get points coordinates.
    x_mr, y_mr, x, y, mag = get_coords(in_file, mag_range)
    # print '\nPhotometric data obtained'
    text += 'Photometric data obtained\n'
    text += 'Total stars: {}\n'.format(len(x))
    text += 'Filtered by {} <= mag < {}: {arg3}\n\n'.format(
        *mag_range, arg3=len(x_mr))

    # Obtain voronoi diagram using the filtered coordinates.
    points = zip(x_mr, y_mr)
    vor = Voronoi(points)
    print 'Voronoi diagram obtained'

    # Get data on Voronoi diagram.
    print '\nProcessing Voronoi diagram'
    acp_pts, rej_pts, pts_area, pts_vert = get_vor_data(points, vor)

    # Obtain average area *using filtered stars*.
    avr_area = ((max(x_mr) - min(x_mr)) * (max(y_mr) - min(y_mr))) / len(x_mr)
    # Filter out large border values.
    pts_area_filt = []
    for _ in pts_area:
        if _ < (5 * avr_area):
            pts_area_filt.append(_)
    # Generate area histogram plot.
    area_hist(f_name, mag_range, pts_area_filt, avr_area)

    # Calculate most probable area for the Voronoi cells.
    most_prob_a = (5. / 7.) * np.mean(pts_area_filt)
    text += ("Most probable area for Voronoi cells (stars filtered by "
             "magnitude): {:.2f}\n".format(most_prob_a))

    # Smaller values imply faster processing since more stars are filtered out.
    for a_f in avr_area_frac:
        pts_thres, vert_thres = area_filter(acp_pts, pts_area, pts_vert,
                                            most_prob_a, a_f)
        text1 = ("\n* Points above ({:.2f} * average_area) threshold: "
                 "{}\n".format(a_f, len(pts_thres)))
        print text1
        text += text1

        print 'Detect shared vertex'
        all_groups = shared_vertex(vert_thres)

        print '\nAssign points to groups'
        neighbors = group_stars(pts_thres, vert_thres, all_groups)

        # Define several minimum neighbors thresholds.
        # for m_n in min_neighbors:
        cent_rad_all = [[], []]
        for m_n in np.arange(300, 9, -10):

            pts_neighbors = neighbors_filter(neighbors, m_n)

            # Check if at least one group was defined with the minimum
            # required number of neighbors.
            if len(pts_neighbors) > 0:
                text += "Groups with more than {} members: {}\n".format(
                    m_n, len(pts_neighbors))

                # Obtain center and radius for each overdensity identified.
                print ("\nObtain center and radius for {} groups with "
                       " {} min neighbors".format(len(pts_neighbors), m_n))
                cent_rad = get_cent_rad(pts_neighbors)

                # Check that the detected overdensities are unique and not
                # found previously.
                new_cent_rads, cent_rad_all[0] = check_cent_rad(
                    cent_rad, cent_rad_all[0])

                # If any *new* overdensities were found, add them to the list.
                if new_cent_rads:
                    cent_rad_all[0] += new_cent_rads
                    cent_rad_all[1] += [m_n] * len(new_cent_rads)

            else:
                text += "No groups with more than {} members.\n".format(
                    m_n)

        if cent_rad_all[0]:
            # Write data to file.
            save_cent_rad(f_name, a_f, cent_rad_all)
            # Plot diagram.
            vor_plot(f_name, a_f, x, y, mag, pts_thres, pts_neighbors,
                     cent_rad_all[0], vor)
        else:
            print 'No groups found with any number of members.'

    # Store info in log file.
    save_to_log(f_name, text, 1)


if __name__ == "__main__":
    start = time.time()
    main()
    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print '\nFull run completed in {:.0f}m {:.0f}s.'.format(m, s)

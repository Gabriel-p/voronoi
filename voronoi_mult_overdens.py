

from scipy.spatial import Voronoi
import time
import numpy as np
from functions.get_params_in import get_params_in
from functions.get_coords import get_coords
from functions.get_cent_rad import get_cent_rad
from functions.group_overdens import merge_overdens
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
    Obtain from Voronoi diagram:
    - acp_pts (points whose diagram/cell/polygon contains no negative indexes)
    - rej_pts (points with negative index vertex)
    - pts_area (area of each cell associated with each accepted point)
    - pts_vert (indexes of vertices associated with each point).
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


def area_filter(acp_pts, pts_area, pts_vert, most_prob_a, a_f):
    '''
    Filter *out* those points whose area is *above* a certain fraction of the
    most probable area value passed.

    http://stackoverflow.com/a/32058576/1391441
    '''

    pts_thres, vert_thres = [], []
    for i, p in enumerate(acp_pts):
        # Keep point if its area is below the maximum threshold.
        if pts_area[i] < most_prob_a * a_f:
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


def group_stars(pts_thres, vert_thres, all_groups):
    '''
    For each defined group, find the points that correspond to it.
    '''

    tot_sols = element_count(all_groups, pts_thres)
    milestones, c = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 0

    neighbors = [[] for _ in range(len(all_groups))]
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

                # Print percentage done.
                milestones = print_perc(c, tot_sols, milestones)

    return neighbors


def neighbors_filter(neighbors, min_neighbors):
    '''
    Filter the points that passed the area filter, discarding those which
    have fewer neighbors than the 'min_neighbors' value.
    '''
    pts_neighbors = []
    for g in neighbors:
        if len(g) > min_neighbors:
            pts_neighbors.append(zip(*g))

    return pts_neighbors


def main():
    # Read parameters from input file.
    in_file, mag_range, avr_area_frac, m_n = get_params_in()

    # Each sub-list in 'in_file' is a row of the file.
    f_name = in_file[:-4]
    # Create log file.
    text = 'Processing file: {}\n'.format(f_name)
    print text
    save_to_log(f_name, text, mag_range[1], 0)
    text = ''

    # Get points coordinates.
    x_mr, y_mr, x, y, mag = get_coords(in_file, mag_range)
    text += 'Photometric data obtained\n'
    text += 'Total stars: {}\n'.format(len(x))
    text1 = 'Stars filtered by {} <= mag < {}: {arg3}\n'.format(
        *mag_range, arg3=len(x_mr))
    print text1
    text += text1

    # Obtain Voronoi diagram using the *magnitude filtered coordinates*.
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

    # Calculate most probable area for the Voronoi cells.
    most_prob_a = (5. / 7.) * np.mean(pts_area_filt)
    text += ("\nMost probable area for Voronoi cells (stars filtered by "
             "magnitude): {:.2f}\n".format(most_prob_a))

    # Smaller values imply faster processing since more stars are filtered out.
    for a_f in avr_area_frac:

        # Generate area histogram plot.
        area_hist(f_name, mag_range, a_f, pts_area_filt, avr_area)

        # Apply maximum area filter.
        pts_thres, vert_thres = area_filter(acp_pts, pts_area, pts_vert,
                                            most_prob_a, a_f)
        text1 = ("\n* Stars below ({:.2f} * most_prob_area) threshold: "
                 "{}\n".format(a_f, len(pts_thres)))
        print text1
        text += text1

        print 'Detect shared vertex'
        all_groups = shared_vertex(vert_thres)

        # This is a costly process.
        print '\nAssign points to groups'
        neighbors = group_stars(pts_thres, vert_thres, all_groups)

        # Keep only those groups with a higher number of members than
        # min_neighbors.
        print '\nDiscard groups with total number of members below the limit.'
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

            old_cent_rad, new_cent_rad = merge_overdens(cent_rad)
            if old_cent_rad:
                text += ('\n{} groups were merged'.format(
                         len(cent_rad) - len(new_cent_rad)))
            else:
                text += '\nNo groups were merged.'
        else:
            new_cent_rad = []
            text += "\nNo groups with more than {} members found.\n".format(
                m_n)

        if new_cent_rad:
            # Write data to file.
            save_cent_rad(f_name, mag_range[1], a_f, m_n, new_cent_rad)
            # Plot diagram.
            print '\nPlotting.'
            vor_plot(f_name, mag_range[1], a_f, m_n, x, y, mag, x_mr, y_mr,
                     pts_thres, pts_neighbors, old_cent_rad, new_cent_rad, vor)
        else:
            print '\nNo groups found with any number of members.'

    # Store info in log file.
    save_to_log(f_name, text, mag_range[1], 1)
    print '\nDone.'


if __name__ == "__main__":
    start = time.time()
    main()
    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print '\nFull run completed in {:.0f}m {:.0f}s.'.format(m, s)

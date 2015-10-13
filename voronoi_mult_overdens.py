

from scipy.spatial import Voronoi, voronoi_plot_2d
import time
import numpy as np
from functions.get_params_in import get_params_in
from functions.get_coords import get_coords
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
    print 'Points above ({} * average area) threshold: {}'.format(
        avr_area_frac, len(pts_thres))

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


def get_cent_rad(pts_neighbors):
    '''
    Assign a center and a radius for each overdensity identified.
    '''
    cent_rad = []
    for g in pts_neighbors:
        cent_x, cent_y = np.mean(g[0]), np.mean(g[1])
        dist_cent = np.sqrt((cent_x - np.array(g[0])) ** 2 +
                            (cent_y - np.array(g[1])) ** 2)
        rad = max(dist_cent)
        cent_rad.append([cent_x, cent_y, rad])

    return cent_rad


def save_to_file(f_name, a_f, m_n, cent_rad):
    '''
    Save center and radius data to file for each parameter value processed.
    '''

    data_out = 'out_fig_dat/' + f_name + '_' + str(a_f) + '_' + \
        str(int(m_n)) + '.out'
    with open(data_out, 'w') as f:
        f.write("#area_frac   min_neigh    cent_x    cent_y    rad\n")
    with open(data_out, 'a') as f:
        for l in cent_rad:
            f.write("{:.3f}{:>11.0f}{:>16.4f}{:>10.4f}{:>10.4f}\n".format(
                a_f, m_n, *l))


def main():
    start = time.time()

    # Read parameters from input file.
    in_file, avr_area_frac, min_neighbors = get_params_in()
    print avr_area_frac
    print min_neighbors
    raw_input()

    # Each sub-list in 'in_file' is a row of the file.
    f_name = in_file[:-4]
    print 'Processing file: {}'.format(f_name)

    # Get points coordinates.
    x, y, mag = get_coords(in_file)
    xmin, xmax, ymin, ymax = min(x), max(x), min(y), max(y)
    avr_area = ((xmax - xmin) * (ymax - ymin)) / len(x)
    print '\nPhotometric data obtained'

    # Obtain voronoi diagram.
    points = zip(x, y)
    vor = Voronoi(points)
    print '\nVoronoi diagram obtained'

    # Get data on Voronoi diagram.
    print '\nProcessing Voronoi diagram'
    acp_pts, rej_pts, pts_area, pts_vert = get_vor_data(points, vor)
    print 'Voronoi diagram processed'
    # Generate area histogram plot.
    area_hist(f_name, pts_area, avr_area)

    # Get overdensities according to the two parameter values.

    # Smaller values imply faster processing since more stars are filtered out.
    for a_f in avr_area_frac:
        pts_thres, vert_thres = area_filter(acp_pts, pts_area, pts_vert,
                                            avr_area, a_f)
        print '\nArea filter applied for {} value'.format(a_f)

        print '\nDetect shared vertex'
        all_groups = shared_vertex(vert_thres)
        print 'Shared vertex processed'

        print '\nAssign points to groups'
        neighbors = group_stars(pts_thres, vert_thres, all_groups)
        print 'Points grouped'

        # Defining several minimum neighbors thresholds is cheap.
        for m_n in min_neighbors:

            pts_neighbors = neighbors_filter(neighbors, m_n)
            print "\nGroups with more than {} members: {}".format(
                m_n, len(pts_neighbors))

            # Obtain center and radius for each overdensity identified.
            cent_rad = get_cent_rad(pts_neighbors)
            print 'Center and radii assigned'

            if len(pts_neighbors) > 0:

                # Write data to file.
                save_to_file(f_name, a_f, m_n, cent_rad)
                print 'Data saved to file'
                # Plot diagram.
                vor_plot(f_name, a_f, m_n, x, y, mag, pts_thres, pts_neighbors,
                         cent_rad)

    # End of run.
    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    print '\nFull iteration completed in {:.0f}m {:.0f}s.'.format(m, s)


if __name__ == "__main__":
    main()

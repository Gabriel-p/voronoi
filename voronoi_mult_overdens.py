

from scipy.spatial import Voronoi
import time
import numpy as np
from functions.percent_done import print_perc
from functions.get_params_in import get_params_in
from functions.get_data import get_data
from functions.get_cent_rad import get_cent_rad
from functions.filt_density import filt_density
from functions.integ_mag import filt_integ_mag
from functions.save_to_file import save_cent_rad, save_to_log
from functions.make_plots import all_plots


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


def get_vor_data(points, mag_mr, vor):
    '''
    Obtain from Voronoi diagram:
    - acp_pts (points whose diagram/cell/polygon contains no negative indexes)
    - rej_pts (points with negative index vertex)
    - pts_area (area of each cell associated with each accepted point)
    - pts_vert (indexes of vertices associated with each point).
    '''

    tot_sols = len(vor.regions)
    milestones = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    acp_pts, acp_mags, rej_pts, pts_area, pts_vert = [], [], [], [], []
    for i, reg in enumerate(vor.regions):
        # If region is not empty.
        if reg:

            # Get index of the point that corresponds to this region.
            p_idx = np.where(vor.point_region == i)[0][0]

            # Check that the vertices index for the region are all positive.
            if all([_ >= 0 for _ in reg]):
                # Store point coordinates.
                acp_pts.append(points[p_idx])
                acp_mags.append(mag_mr[p_idx])

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

    return acp_pts, acp_mags, rej_pts, pts_area, pts_vert


def large_area_filt(pts_area, avr_area):
    '''
    Filter out polygons located at the borders of the frame with very large
    area values.
    Find mean polygon area for the Voronoi cells of the Poisson distribution.
    '''
    pts_area_filt = []
    for _ in pts_area:
        if _ < (5 * avr_area):
            pts_area_filt.append(_)

    # Calculate mean area for all the Voronoi cells.
    mean_filt_area = np.mean(pts_area_filt)

    return pts_area_filt, mean_filt_area


def main():
    start = time.time()

    # Read parameters from input file.
    in_file, in_file_cols, coords_flag, vor_flag, cr_file_cols, mag_range,\
        area_frac_range, m_m_n, intens_frac, dens_frac = get_params_in()

    # Each sub-list in 'in_file' is a row of the file.
    f_name = in_file[:-4]
    # Create log file.
    save_to_log(f_name, 'Processing file: {}'.format(f_name), 'w')
    save_to_log(f_name, 'Mag range: {}'.format(mag_range), 'a')
    if vor_flag == 'voronoi':
        save_to_log(f_name, 'Area range: {}'.format(area_frac_range), 'a')
    save_to_log(f_name, 'Min/max neighbors: {}'.format(m_m_n), 'a')
    save_to_log(f_name, "Frame's density fraction: {}".format(
        intens_frac), 'a')
    save_to_log(f_name, "Frame's I/A fraction: {}\n".format(dens_frac), 'a')

    # Get points coordinates and magnitudes.
    x_mr, y_mr, mag_mr, x, y, mag, ra_cent, dec_cent = get_data(
        in_file, in_file_cols, mag_range, coords_flag, vor_flag)
    save_to_log(f_name, 'Photometric data obtained', 'a')
    save_to_log(f_name, 'Total number of stars: {}'.format(len(x)), 'a')
    if coords_flag == 'deg':
        save_to_log(f_name, 'Center of frame: RA={}, DEC={}'.format(
            ra_cent, dec_cent), 'a')
    else:
        save_to_log(f_name, 'Center of frame: x={}, y={}'.format(
            ra_cent, dec_cent), 'a')
    save_to_log(f_name, '\nStars filtered by mag range: {}'.format(len(x_mr)),
                'a')

    # Obtain Voronoi diagram using the *magnitude filtered coordinates*.
    points = zip(x_mr, y_mr)
    vor = Voronoi(points)
    save_to_log(f_name, 'Voronoi diagram obtained', 'a')

    # Get data on Voronoi diagram.
    save_to_log(f_name, 'Processing Voronoi diagram', 'a')
    acp_pts, acp_mags, rej_pts, pts_area, pts_vert = get_vor_data(points,
                                                                  mag_mr, vor)

    # Obtain average area *using magnitude filtered stars*.
    fr_area = ((max(x_mr) - min(x_mr)) * (max(y_mr) - min(y_mr)))
    avr_area = fr_area / len(x_mr)

    # Filter out large area values for border polygons.
    pts_area_filt, mean_filt_area = large_area_filt(pts_area, avr_area)
    save_to_log(f_name, ("Mean Voronoi cells area (stars in mag range): "
                         "{:g} {}^2".format(mean_filt_area, coords_flag)), 'a')

    # Obtain center coordinates and radii either automatically by grouping
    # neighbor stars, or by reading those values from a file.
    # This function applies the area and neighbors filters.
    cent_rad, pts_area_thres, mag_area_thres, pts_neighbors,\
        mags_neighbors = get_cent_rad(
            f_name, coords_flag, cr_file_cols, m_m_n, vor_flag, ra_cent,
            dec_cent, acp_pts, acp_mags, pts_area, pts_vert, mean_filt_area,
            area_frac_range)

    # Filter/organize groups found by their density, using stars filtered
    # by magnitude.
    dens_accp_groups, dens_rej_groups, area_frac_range = filt_density(
        f_name, fr_area, x_mr, y_mr, mag_mr, cent_rad, vor_flag,
        area_frac_range, dens_frac, mean_filt_area)
    save_to_log(
        f_name, "\nGroups filtered by density (stars/area): {}".format(
            len(dens_accp_groups)), 'a')

    # Filter/organize groups found by their I/A, using stars filtered
    # by magnitude.
    intens_acc_dens_acc, intens_acc_dens_rej, intens_rej_dens_acc,\
        intens_rej_dens_rej = filt_integ_mag(
            x_mr, y_mr, mag_mr, dens_accp_groups, dens_rej_groups,
            intens_frac)
    save_to_log(
        f_name, "\nGroups filtered by intensity/area: {}".format(
            len(intens_acc_dens_acc[0])), 'a')

    # Write data to file.
    save_cent_rad(f_name, cent_rad, intens_acc_dens_acc,
                  intens_acc_dens_rej, intens_rej_dens_acc,
                  intens_rej_dens_rej)
    save_to_log(f_name, "\nData saved to {}.out file".format(f_name), 'a')

    # Make plots.
    save_to_log(f_name, '\nPlotting', 'a')
    all_plots(f_name, mag_range, area_frac_range, x, y, mag,
              coords_flag, vor_flag, x_mr, y_mr, mag_mr, pts_area_filt,
              pts_area_thres, mag_area_thres, pts_neighbors, mags_neighbors,
              intens_frac, dens_frac, cent_rad, intens_acc_dens_acc,
              intens_acc_dens_rej, intens_rej_dens_acc, intens_rej_dens_rej)

    # Done.
    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    save_to_log(f_name, 'Full run completed in {:.0f}m {:.0f}s.'.format(
        m, s), 'a')


if __name__ == "__main__":
    main()

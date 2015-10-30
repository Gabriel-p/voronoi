
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.offsetbox as offsetbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage.filters import gaussian_filter
import bisect
# from scipy.spatial import voronoi_plot_2d


def save_plot(f_name, fig_id, fig):
    '''
    Save output png file.
    '''
    fig.tight_layout()
    fig_name = 'out_fig_dat/' + f_name + '_' + fig_id + '.png'
    plt.savefig(fig_name, dpi=300)


def vor_2d_poisson(y):
    '''
    Function that fits the normalized area of Voronoi cells for a 2D Poisson
    distribution.
    Taken from: On the size-distribution of Poisson Voronoi cells,
    Jarai-Szabo & Neda (2008); http://arxiv.org/abs/cond-mat/0406116 (Eq. 10)
    '''
    # y = vor_cell_area / <vor_cell_area>
    a = (343. / 15.) * np.sqrt(7. / (2. * np.pi))
    f_y = a * (y ** (5. / 2.)) * np.exp(-1. * (7. / 2.) * y)

    return f_y


def area_hist(f_name, mag_range, area_frac_range, pts_area_filt):
    '''
    '''
    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(111)
    plt.xlabel('(polygon area) / <polygon area>', fontsize=12)
    plt.ylabel("Normalized distribution", fontsize=12)
    # Area of each point's polygon as a fraction of the average area.
    pts_area_f_mean = np.mean(pts_area_filt)
    frac_area = pts_area_filt / pts_area_f_mean
    # Vertical lines.
    plt.axvline(x=(5. / 7.), color='k', ls='--', lw=2,
                label='Most prob polygon area ({:.2E} u^2)'.format(
                    (5. / 7.) * pts_area_f_mean))
    plt.axvline(x=1., color='k', ls='--', lw=2,
                label='Mean Voronoi cell area ({:.2E} u^2)'.format(
                    pts_area_f_mean))
    a_f_min, a_f_max = [round(_, 2) for _ in area_frac_range]
    plt.axvline(x=a_f_min, color='r', ls='--', lw=2,
                label='Min area fraction: {:.2f} ({:.2E} u^2)'.format(
                    a_f_min, a_f_min * pts_area_f_mean))
    plt.axvline(x=a_f_max, color='r', ls='--', lw=2,
                label='Max area fraction: {:.2f} ({:.2E} u^2)'.format(
                    a_f_max, a_f_max * pts_area_f_mean))
    # Normalized histogram.
    weights = np.ones_like(frac_area)/len(frac_area)
    plt.hist(frac_area, color='#C6D8E5', bins=50, range=[0., 4.],
             weights=weights, normed=True)
    # Plot theoretical fit to 2D Poisson Voronoi area distribution.
    x = np.arange(0., 4., 0.05)
    plt.plot(x, vor_2d_poisson(x), c='#117050', lw=3.5,
             label='Theoretical distribution')
    # Add text box.
    text = '{} <= mag < {}'.format(*mag_range)
    ob = offsetbox.AnchoredText(text, pad=0.5, loc=7, prop=dict(size=11))
    ob.patch.set(alpha=0.5)
    ax1.add_artist(ob)
    # get handles
    handles, labels = ax1.get_legend_handles_labels()
    # use them in the legend
    ax1.legend(handles, labels, loc='upper right', numpoints=2, fontsize=11)
    # Save plot to file.
    save_plot(f_name, 'area_histo', fig)


def dens_vs_intens_plot(f_name, mag_range, dens_frac, intens_frac,
                        intens_acc_dens_acc, intens_acc_dens_rej,
                        intens_rej_dens_acc, intens_rej_dens_rej):
    '''
    '''
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111)
    # Clusters that passed both the density and I/A filters.
    x_aa, y_aa = intens_acc_dens_acc[3], intens_acc_dens_acc[2]
    # Clusters that passed the I/A but not the density filter.
    x_ar, y_ar = intens_acc_dens_rej[3], intens_acc_dens_rej[2]
    # Clusters that passed the density but not the I/A filter.
    x_ra, y_ra = intens_rej_dens_acc[3], intens_rej_dens_acc[2]
    # Clusters that passed neither the density nor the I/A filters.
    x_rr, y_rr = intens_rej_dens_rej[3], intens_rej_dens_rej[2]
    max_x = max(x_aa + x_ar + x_ra + x_rr)
    max_y = max(y_aa + y_ar + y_ra + y_rr)
    # Limit axis extension to some reasonable value.
    if max_x > np.mean(x_aa + x_ar + x_ra + x_rr) * 10.:
        max_x = np.mean(x_aa + x_ar + x_ra + x_rr) * 10.
    if max_y > np.mean(y_aa + y_ar + y_ra + y_rr) * 10.:
        max_y = np.mean(y_aa + y_ar + y_ra + y_rr) * 10.
    plt.xlim(0., max_x + max_x * 0.1)
    plt.ylim(0., max_y + max_y * 0.1)
    plt.xlabel("Intensity / area", fontsize=12)
    plt.ylabel("Density (N in mag range) / area", fontsize=12)
    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    plt.grid(b=True, which='minor', color='gray', linestyle='-', zorder=1)
    plt.axhline(y=1., color='goldenrod', ls='--', lw=2,
                label='(N in frame) / area')
    plt.axvline(x=1., color='b', ls='--', lw=2,
                label="(Frame's intensity) / area")
    # Scatter points.
    plt.scatter(x_aa, y_aa, c='g', marker='o', lw=0.2, s=50,
                label='Dens + I/A accp ({})'.format(len(x_aa)), zorder=4)
    plt.scatter(x_ar, y_ar, c='c', marker='o', lw=0.2, s=50,
                label='I/A accp, Dens rej ({})'.format(len(x_ar)), zorder=4)
    plt.scatter(x_ra, y_ra, c='m', marker='o', lw=0.2, s=50,
                label='Dens accp, I/A rej ({})'.format(len(x_ra)), zorder=4)
    plt.scatter(x_rr, y_rr, c='r', marker='o', lw=0.2, s=50,
                label='Dens + I/A rej ({})'.format(len(x_rr)), zorder=4)
    # Add text box.
    text = '     Mag range:\n{} <= mag < {}'.format(*mag_range)
    ob = offsetbox.AnchoredText(text, pad=0.5, loc=1, prop=dict(size=10))
    ob.patch.set(alpha=0.5)
    ax.add_artist(ob)
    # Legend.
    leg = plt.legend(fancybox=True, loc='lower right', scatterpoints=1,
                     fontsize=10, markerscale=1.)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.85)
    # Save plot to file.
    save_plot(f_name, 'dens_vs_intens', fig)


def intensity_plot(f_name, mag_range, intens_frac, intens_acc_dens_acc,
                   intens_rej_dens_acc):
    '''
    '''
    fig = plt.figure(figsize=(20, 10))
    ax1 = plt.subplot(121)
    plt.xlabel("Intensity / area", fontsize=12)
    plt.ylabel("Normalized distribution", fontsize=12)
    # Vertical lines.
    plt.axvline(x=1., color='k', ls='--', lw=2,
                label="(Frame's intensity) / area = 1.")
    plt.axvline(x=intens_frac, color='r', ls='--', lw=2,
                label='(Minimum intensity) / area = {}'.format(intens_frac))
    # Normalized histogram.
    intens_area = intens_acc_dens_acc[3] + intens_rej_dens_acc[3]
    max_val = max(intens_area)
    weights = np.ones_like(intens_area)/len(intens_area)
    plt.hist(intens_area, color='#C6D8E5', bins=50, range=[0., max_val],
             weights=weights, normed=True)
    # Add text box.
    text = '     Mag range:\n{} <= mag < {}'.format(*mag_range)
    ob = offsetbox.AnchoredText(text, pad=0.5, loc=7, prop=dict(size=10))
    ob.patch.set(alpha=0.5)
    ax1.add_artist(ob)
    # get handles
    handles, labels = ax1.get_legend_handles_labels()
    # use them in the legend
    ax1.legend(handles, labels, loc='upper right', numpoints=2, fontsize=11)

    ax2 = plt.subplot(122)
    # Accepted clusters.
    if intens_acc_dens_acc[0]:
        x_a, y_a, z_a = zip(*intens_acc_dens_acc[0])[2],\
            intens_acc_dens_acc[1], intens_acc_dens_acc[3]
        # Order lists to put min rad values on top.
        ord_za, ord_xa, ord_ya = map(list, zip(*sorted(zip(z_a, x_a, y_a),
                                     reverse=False)))
    else:
        ord_za, ord_xa, ord_ya = [], [], []
    # Rejected clusters.
    if intens_rej_dens_acc[0]:
        x_r, y_r, z_r = map(list, (zip(*intens_rej_dens_acc[0])[2],
                            intens_rej_dens_acc[1], intens_rej_dens_acc[3]))
    else:
        x_r, y_r, z_r = [], [], []
    # Limits.
    try:
        max_x, max_y = max(ord_xa + x_r), max(ord_ya + y_r)
        v_max = max(ord_za + z_r)
    except:
        # If all lists are empty, set values to 0.
        max_x, max_y, v_max = 0., 0., 0.
    v_min = 0.
    # Plot.
    plt.xlim(0., max_x + max_x * 0.1)
    plt.ylim(-1, max_y + max_y * 0.1)
    plt.xlabel("Radius", fontsize=12)
    plt.ylabel("N stars (in mag range)", fontsize=12)
    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='gray', linestyle='-', zorder=1)
    plt.grid(b=True, which='minor', color='gray', linestyle='-', zorder=1)
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Only plot if list is not empty.
    if x_r:
        plt.scatter(x_r, y_r, c=z_r, marker='s', lw=0.2, s=35, cmap=cm,
                    vmin=v_min, vmax=v_max,
                    label='Rejected overdensities ({})'.format(len(x_r)),
                    zorder=3)
    SC = plt.scatter(ord_xa, ord_ya, c=ord_za, marker='o', lw=0.2, s=50,
                     cmap=cm, vmin=v_min, vmax=v_max,
                     label='Accepted overdensities ({})'.format(len(ord_xa)),
                     zorder=4)
    # Add text box.
    text = '(Minimum intensity) / area: {}'.format(intens_frac)
    ob = offsetbox.AnchoredText(text, pad=0.5, loc=6, prop=dict(size=10))
    ob.patch.set(alpha=0.5)
    ax2.add_artist(ob)
    # Legend.
    leg = plt.legend(fancybox=True, loc='upper left', scatterpoints=1,
                     fontsize=10, markerscale=1.)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.85)
    # Position colorbar.
    the_divider = make_axes_locatable(ax2)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    cbar.set_label("Intensity / area", fontsize=12, labelpad=5)

    # Save plot to file.
    save_plot(f_name, 'intensity', fig)


def star_size(mag_data, min_mag):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 50.
    return 0.75 + factor * 10 ** ((np.array(mag_data) - min_mag) / -2.5)


def dens_map(x_data, y_data, new_cent_rad):
    '''
    Create a 2D density map according to parameters passed.
    '''
    xmin, xmax = min(x_data), max(x_data)
    ymin, ymax = min(y_data), max(y_data)
    # Calculate the number of bins used.
    x_rang, y_rang = (xmax - xmin), (ymax - ymin)
    bin_width = min(x_rang, y_rang) / 100.
    # Number of bins in x,y given the bin width.
    binsxy = [int(x_rang / bin_width), int(y_rang / bin_width)]

    st_dev = 1.5
    hist, xedges, yedges = np.histogram2d(
        x_data, y_data, range=[[xmin, xmax], [ymin, ymax]], bins=binsxy)
    h_g = gaussian_filter(hist, st_dev, mode='constant')

    dens_cent_rad = []
    for c_r in new_cent_rad:
        x_cent_bin = bisect.bisect_left(xedges, c_r[0])
        y_cent_bin = bisect.bisect_left(yedges, c_r[1])
        # Store center bin coords for the filtered hist.
        dens_cent_rad.append([(x_cent_bin - 1), (y_cent_bin - 1),
                             c_r[2] / bin_width])

    return h_g, dens_cent_rad


def vor_plot(f_name, x, y, mag, coords_flag, vor_flag, x_mr, y_mr, mag_mr,
             pts_area_thres, mag_area_thres, pts_neighbors, mags_neighbors,
             cent_rad, intens_acc_dens_acc,
             intens_acc_dens_rej, intens_rej_dens_acc, intens_rej_dens_rej):
    # figure size.
    fig = plt.figure(figsize=(40, 20))
    gs = gridspec.GridSpec(2, 4)

    ax = plt.subplot(gs[0:2, 0:2])
    # plt.gca().set_aspect('equal')
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    if coords_flag == 'deg':
        # If RA is used, invert axis.
        ax.invert_xaxis()
        x_lab, y_lab = 'RA', 'DEC'
    else:
        x_lab, y_lab = 'X', 'Y'
    plt.xlabel('{} ({})'.format(x_lab, coords_flag), fontsize=12)
    plt.ylabel('{} ({})'.format(y_lab, coords_flag), fontsize=12)
    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='k', linestyle='-', zorder=1)
    plt.grid(b=True, which='minor', color='k', linestyle='-', zorder=1)

    # All points.
    min_mag = min(mag)
    st_sizes_arr = star_size(mag, min_mag)
    plt.scatter(x, y, c='#4C4C4C', marker='o', s=st_sizes_arr, lw=0.,
                label='All stars')
    # plt.scatter(x, y, c='gray', marker='.', s=7, lw=0., label='All stars')
    # Points that passed the magnitude filter.
    st_sizes_arr = star_size(mag_mr, min_mag)
    plt.scatter(x_mr, y_mr, marker='o', s=st_sizes_arr, lw=0.2,
                facecolor='none', edgecolor='b', label='Mag filter')
    # Only plot if these points were obtained.
    if vor_flag == 'voronoi':
        # Points that passed the area threshold filter.
        x_t, y_t = zip(*pts_area_thres)
        st_sizes_arr = star_size(mag_area_thres, min_mag)
        plt.scatter(x_t, y_t, marker='o', s=st_sizes_arr, lw=0.3,
                    facecolor='none', edgecolor='gold', label='Area filter')
        # Neighbor points.
        x_n, y_n = [], []
        for g in pts_neighbors:
            x_n += list(g[0])
            y_n += list(g[1])
        st_sizes_arr = star_size(mags_neighbors, min_mag)
        plt.scatter(x_n, y_n, marker='o', s=st_sizes_arr, lw=0.4,
                    facecolor='none', edgecolor='g', label='Neighbor stars')
    # Legend.
    leg = plt.legend(fancybox=True, loc='upper right', scatterpoints=1,
                     fontsize=10, markerscale=1.)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.85)

    # Print radii for all detected groups above the minimum neighbors limit.
    for c_r in cent_rad:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='b', ls='dashed',
                            fill=False, lw=1.2)
        fig.gca().add_artist(circle)

    # Print radii for dens+I/A accepted groups.
    for c_r in intens_acc_dens_acc[0]:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='g',
                            fill=False, lw=1.2)
        fig.gca().add_artist(circle)
    # Print radii for dens acc & I/A rej groups.
    for c_r in intens_rej_dens_acc[0]:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='m',
                            fill=False, lw=1.2)
        fig.gca().add_artist(circle)
    # Print radii for dens rej & I/A acc groups.
    for c_r in intens_acc_dens_rej[0]:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='c',
                            fill=False, lw=1.2)
        fig.gca().add_artist(circle)
    # Print radii for dens+I/A rejected groups.
    for c_r in intens_rej_dens_rej[0]:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='r',
                            fill=False, lw=1.2)
        fig.gca().add_artist(circle)

    # # Plot Voronoi cells
    # # voronoi_plot_2d(vor)
    # plt.plot(vor.vertices[:, 0], vor.vertices[:, 1], '+')
    # for simplex in vor.ridge_vertices:
    #     simplex = np.asarray(simplex)
    #     if np.all(simplex >= 0):
    #         plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1],
    #              'k-', lw=0.5)

    # x, y = zip(*rej_pts)
    # plt.scatter(x, y, s=70, c='r', zorder=3)

    # x, y = zip(*acp_pts)
    # plt.scatter(x, y, s=50, c='g', marker='s')

    # x, y = zip(*pts_thres)
    # plt.scatter(x, y, s=30, c='b', marker='o')

    cm = plt.cm.gist_earth_r

    # All stars density map.
    ax1 = plt.subplot(gs[0, 2])
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    h_g, dens_cent_rad = dens_map(x, y, intens_acc_dens_acc[0])
    ax1.imshow(h_g.transpose(), origin='lower', cmap=cm, aspect='auto')
    # Print radii
    for c_r in dens_cent_rad:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='r', lw=1.5,
                            fill=False)
        fig.gca().add_artist(circle)
    # Add text box.
    text = '1- All stars'
    ob = offsetbox.AnchoredText(text, pad=0.5, loc=1, prop=dict(size=10))
    ob.patch.set(alpha=0.5)
    ax1.add_artist(ob)

    # Magnitude filtered stars density map.
    ax2 = plt.subplot(gs[0, 3])
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)
    # Add extreme points so aspect is the same for all density maps.
    x_mr = x_mr + [min(x), max(x)]
    y_mr = y_mr + [min(y), max(y)]
    h_g, dens_cent_rad = dens_map(x_mr, y_mr, intens_acc_dens_acc[0])
    ax2.imshow(h_g.transpose(), origin='lower', cmap=cm, aspect='auto')
    # Print radii
    for c_r in dens_cent_rad:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='r', lw=1.5,
                            fill=False)
        fig.gca().add_artist(circle)
    # Add text box.
    text = '2- Magnitude filter'
    ob = offsetbox.AnchoredText(text, pad=0.5, loc=1, prop=dict(size=10))
    ob.patch.set(alpha=0.5)
    ax2.add_artist(ob)

    # Only plot if these points were obtained.
    if vor_flag == 'voronoi':
        # Area threshold filtered stars density map.
        ax3 = plt.subplot(gs[1, 2])
        ax3.xaxis.set_visible(False)
        ax3.yaxis.set_visible(False)
        # Add extreme points so aspect is the same for all density maps.
        x_t = list(x_t) + [min(x), max(x)]
        y_t = list(y_t) + [min(y), max(y)]
        h_g, dens_cent_rad = dens_map(x_t, y_t, intens_acc_dens_acc[0])
        ax3.imshow(h_g.transpose(), origin='lower', cmap=cm, aspect='auto')
        # Print radii
        for c_r in dens_cent_rad:
            circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='r', lw=1.5,
                                fill=False)
            fig.gca().add_artist(circle)
        # Add text box.
        text = '3- Area filter'
        ob = offsetbox.AnchoredText(text, pad=0.5, loc=1, prop=dict(size=10))
        ob.patch.set(alpha=0.5)
        ax3.add_artist(ob)

        # Neighbors stars density map.
        ax4 = plt.subplot(gs[1, 3])
        # ax4.set_aspect(aspect=1)
        ax4.xaxis.set_visible(False)
        ax4.yaxis.set_visible(False)
        # Add extreme points so aspect is the same for all density maps.
        x_n = list(x_n) + [min(x), max(x)]
        y_n = list(y_n) + [min(y), max(y)]
        h_g, dens_cent_rad = dens_map(x_n, y_n, intens_acc_dens_acc[0])
        ax4.imshow(h_g.transpose(), origin='lower', cmap=cm, aspect='auto')
        # Print radii
        for c_r in dens_cent_rad:
            circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='r', lw=1.5,
                                fill=False)
            fig.gca().add_artist(circle)
        # Add text box.
        text = '4- Neighbor stars'
        ob = offsetbox.AnchoredText(text, pad=0.5, loc=1, prop=dict(size=10))
        ob.patch.set(alpha=0.5)
        ax4.add_artist(ob)

    # Save plot to file.
    save_plot(f_name, 'voronoi', fig)


def all_plots(f_name, mag_range, area_frac_range, x, y, mag,
              coords_flag, vor_flag, x_mr, y_mr, mag_mr, pts_area_filt,
              pts_area_thres, mag_area_thres, pts_neighbors, mags_neighbors,
              intens_frac, dens_frac, cent_rad, intens_acc_dens_acc,
              intens_acc_dens_rej, intens_rej_dens_acc, intens_rej_dens_rej):
    '''
    Make all plots.
    '''
    # Generate area histogram plot.
    area_hist(f_name, mag_range, area_frac_range, pts_area_filt)

    if pts_neighbors:
        # Density versus intensities plots.
        dens_vs_intens_plot(f_name, mag_range, dens_frac, intens_frac,
                            intens_acc_dens_acc, intens_acc_dens_rej,
                            intens_rej_dens_acc, intens_rej_dens_rej)

        # Intensities plots.
        intensity_plot(f_name, mag_range, intens_frac, intens_acc_dens_acc,
                       intens_rej_dens_acc)

        # Plot main diagram.
        vor_plot(f_name, x, y, mag, coords_flag, vor_flag, x_mr, y_mr, mag_mr,
                 pts_area_thres, mag_area_thres, pts_neighbors, mags_neighbors,
                 cent_rad, intens_acc_dens_acc,
                 intens_acc_dens_rej, intens_rej_dens_acc, intens_rej_dens_rej)

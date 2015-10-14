
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from scipy.spatial import voronoi_plot_2d


def save_plot(f_name, fig_id, fig, a_f):
    '''
    Save output png file.
    '''
    fig.tight_layout()
    if a_f == '':
        fig_name = 'out_fig_dat/' + f_name + '_' + fig_id + '.png'
    else:
        fig_name = 'out_fig_dat/' + f_name + '_' + fig_id + '_' + str(a_f) +\
            '_' + '.png'
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


def area_hist(f_name, mag_range, pts_area_filt, avr_area):
    '''
    '''
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111)
    plt.xlabel('(polygon area) / <polygon area>', fontsize=12)
    plt.ylabel("Voronoi cell areas (normalized)", fontsize=12)
    # Area of each point's polygon as a fraction of the average area.
    pts_area_f_mean = np.mean(pts_area_filt)
    frac_area = pts_area_filt / pts_area_f_mean
    # Normalized histogram.
    weights = np.ones_like(frac_area)/len(frac_area)
    plt.hist(frac_area, color='#C6D8E5', bins=50, range=[0., 4.],
             weights=weights, normed=True)
    # Plot theoretical fit to 2D Poisson Voronoi area distribution.
    x = np.arange(0., 4., 0.05)
    plt.plot(x, vor_2d_poisson(x), c='#117050', lw=3.5)
    plt.axvline(x=(avr_area / pts_area_f_mean), color='r', ls='--',
                lw=2,
                label='<star area> / <polygon area>\n{:.2f} / {:.2f}'.format(
                    avr_area, pts_area_f_mean))
    most_prob_a = (5. / 7.) * pts_area_f_mean
    plt.axvline(x=(5. / 7.), color='k', ls='--', lw=2,
                label='Most probable\npolygon area ({:.2f})'.format(
                    most_prob_a))
    # Add text box.
    text = '{} <= mag < {}'.format(*mag_range)
    ob = offsetbox.AnchoredText(text, pad=0.5, loc=7, prop=dict(size=13))
    ob.patch.set(alpha=0.5)
    ax.add_artist(ob)
    # get handles
    handles, labels = ax.get_legend_handles_labels()
    # use them in the legend
    ax.legend(handles, labels, loc='upper right', numpoints=2, fontsize=13)
    # Save plot to file.
    save_plot(f_name, 'area_histo', fig, '')


def star_size(mag_data):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 500. * (1 - 1 / (1 + 150 / len(mag_data) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag_data) - min(mag_data)) / -2.5)


def vor_plot(f_name, a_f, x, y, mag, pts_thres, neighbors_plot, cent_rad, vor):
    # figure size.
    fig = plt.figure(figsize=(20, 20))

    plt.gca().set_aspect('equal')
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    # plt.gca().invert_xaxis()
    # plt.xlabel('{} ({})'.format('RA', 'deg'), fontsize=12)
    # plt.ylabel('{} ({})'.format('DEC', 'deg'), fontsize=12)
    plt.xlabel('{} ({})'.format('X', 'px'), fontsize=12)
    plt.ylabel('{} ({})'.format('Y', 'px'), fontsize=12)

    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='k', linestyle='-', zorder=1)
    plt.grid(b=True, which='minor', color='k', linestyle='-', zorder=1)
    # st_sizes_arr = star_size(mag)
    # plt.scatter(x, y, c='k', s=st_sizes_arr)
    plt.scatter(x, y, c='k', marker='.', s=1)

    # voronoi_plot_2d(vor)
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

    for g in neighbors_plot:
        plt.scatter(g[0], g[1], s=10, c='c', marker='x')

    # Print radii
    for c_r in cent_rad:
        circle = plt.Circle((c_r[0], c_r[1]), c_r[2], color='r',
                            fill=False, lw=1.5)
        fig.gca().add_artist(circle)

    # Save plot to file.
    save_plot(f_name, 'voronoi', fig, a_f)

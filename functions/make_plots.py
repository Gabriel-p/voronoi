
import numpy as np
import matplotlib.pyplot as plt


def save_plot(f_name, fig_id, fig, a_f, m_n):
    '''
    Save output png file.
    '''
    fig.tight_layout()
    if a_f == '' and m_n == '':
        fig_name = 'out_fig_dat/' + f_name + '_' + fig_id + '.png'
    else:
        fig_name = 'out_fig_dat/' + f_name + '_' + fig_id + '_' + str(a_f) +\
            '_' + str(int(m_n)) + '.png'
    plt.savefig(fig_name, dpi=300)


def area_hist(f_name, pts_area, avr_area):
    '''
    '''
    fig = plt.figure(figsize=(10, 10))
    plt.hist(pts_area / avr_area, bins=50, range=[0., 3.])
    plt.axvline(x=1., color='r', ls='-')

    # Save plot to file.
    save_plot(f_name, 'area_histo', fig, '', '')


def star_size(mag_data):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    factor = 500. * (1 - 1 / (1 + 150 / len(mag_data) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag_data) - min(mag_data)) / -2.5)


def vor_plot(f_name, a_f, m_n, x, y, mag, pts_thres, neighbors_plot, cent_rad):
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
    save_plot(f_name, 'voronoi', fig, a_f, m_n)

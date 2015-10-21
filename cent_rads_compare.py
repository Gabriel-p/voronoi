
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def save_plot(f_name, fig_id, fig):
    '''
    Save output png file.
    '''
    fig.tight_layout()
    fig_name = 'out_fig_dat/' + f_name + '_' + fig_id + '.png'
    plt.savefig(fig_name, dpi=300)


def make_plot(data1, data2):
    '''
    '''
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(111)
    x_min, x_max = 0., 6400
    y_min, y_max = 0., 6400
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.xlabel('{} ({})'.format('RA', 'px'), fontsize=12)
    plt.ylabel('{} ({})'.format('DEC', 'px'), fontsize=12)
    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5,
             zorder=1)
    plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.5,
             zorder=1)
    # Plot Nort-East arrows
    plt.arrow(50, 160, 0, -100, head_width=40, head_length=25, fc='k', ec='k', zorder=3)
    plt.arrow(50, 160, 100, 0, head_width=40, head_length=25, fc='k', ec='k', zorder=3)
    # Print radii for Voronoi groups.
    for c_x, c_y, r in data1:
        circle = plt.Circle((c_x, c_y), r, color='g', fill=False, lw=0.5,
                            zorder=3)
        fig.gca().add_artist(circle)

    # Print radii for final groups.
    # for cc_x, cc_y, cc_r in zip(*[cc_x, cc_y, cc_r]):
    for cc_x, cc_y, cr_x, cr_y in data2:

        # Plot ellipsoid.
        ellipse = Ellipse(xy=(cc_x, cc_y), width=cr_x, height=cr_y,
                          edgecolor='b', angle=0., fc='None', lw=0.5, zorder=3)
        fig.gca().add_artist(ellipse)

        # circle = plt.Circle((cc_x, cc_y), cc_r, color='b', ls='dashed',
        #                     fill=False, lw=0.5, zorder=3)
        # fig.gca().add_artist(circle)

    # Save plot to file.
    save_plot('junk', 'bica', fig)


# Load Voronoi centers and radii.
vor_file = 'out_fig_dat/junk_17.5_0.45_1.0_5.out'
file_data = np.loadtxt(vor_file)
i, j, k = 4, 3, 5
# Extract coordinates and zip them into lists.
c_x, c_y, r = zip(*file_data)[i], zip(*file_data)[j],  zip(*file_data)[k]

# Load compare file centers and radii.
com_file = 'bica.dat'
file_data = np.loadtxt(com_file)
i, j, k, q = 0, 1, 2, 3
# Extract coordinates and zip them into lists.
cc_x, cc_y, cr_x, cr_y = zip(*file_data)[i], zip(*file_data)[j], \
    zip(*file_data)[k], zip(*file_data)[q]
cc_r = (np.array(cr_x) + np.array(cr_y)) / 2.

# Make plot.
data1, data2 = zip(*[c_x, c_y, r]), zip(*[cc_x, cc_y, cr_x, cr_y])
make_plot(data1, data2)

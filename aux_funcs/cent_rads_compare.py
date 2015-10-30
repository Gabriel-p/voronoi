
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def save_plot(f_name, fig_id, fig):
    '''
    Save output png file.
    '''
    fig.tight_layout()
    fig_name = '../out_fig_dat/' + f_name + '_' + fig_id + '.png'
    plt.savefig(fig_name, dpi=300)


def star_size(mag_data):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    # factor = 2000. * (1 - 1 / (1 + 150 / len(mag_data) ** 0.85))
    factor = 50.
    return 0.3 + factor * 10 ** ((np.array(mag_data) - min(mag_data)) / -2.5)


def make_plot(x, y, mag, data1, data2):
    '''
    '''
    fig = plt.figure(figsize=(20, 20))
    ax = plt.subplot(111)
    # plt.gca().set_aspect('equal')
    # x_min, x_max = 0., 6400
    # y_min, y_max = 0., 6400
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # If RA is used, invert axis.
    ax.invert_xaxis()
    # ax.invert_yaxis()
    plt.xlabel('{} ({})'.format('RA', 'deg'), fontsize=12)
    plt.ylabel('{} ({})'.format('DEC', 'deg'), fontsize=12)
    # Set minor ticks
    plt.minorticks_on()
    # Set grid
    plt.grid(b=True, which='major', color='gray', linestyle='-', lw=0.5,
             zorder=1)
    plt.grid(b=True, which='minor', color='gray', linestyle='-', lw=0.5,
             zorder=1)
    # Plot Nort-East arrows
    plt.arrow(50, 160, 0, -100, head_width=40, head_length=25, fc='k', ec='k',
              zorder=3)
    plt.arrow(50, 160, 100, 0, head_width=40, head_length=25, fc='k', ec='k',
              zorder=3)
    # Plot all stars.
    st_sizes_arr = star_size(mag)
    plt.scatter(x, y, c='k', s=st_sizes_arr, lw=0., marker='o')
    # Print radii for Voronoi groups.
    for c_x, c_y, r in data1:
        circle = plt.Circle((c_x, c_y), r, color='g', fill=False, lw=0.5,
                            zorder=3)
        fig.gca().add_artist(circle)
    # Print radii for final groups.
    # for cc_x, cc_y, cc_r in zip(*[cc_x, cc_y, cc_r]):
    for cc_x, cc_y, cr_x, cr_y in data2:

        # Plot ellipsoid. Multiply by 2 because 'width' and 'height' are
        # total extensions.
        ellipse = Ellipse(xy=(cc_x, cc_y), width=cr_x * 2, height=cr_y * 2,
                          edgecolor='r', angle=0., fc='None', lw=0.5, zorder=3)
        fig.gca().add_artist(ellipse)

        # circle = plt.Circle((cc_x, cc_y), cc_r, color='b', ls='dashed',
        #                     fill=False, lw=0.5, zorder=3)
        # fig.gca().add_artist(circle)

    # Save plot to file.
    save_plot('junk', 'bica', fig)


# Load all data.
dat_file = '../junk_psf.dat'
file_data = np.loadtxt(dat_file)
i, j, k = 0, 1, 8
print 'All stars: {}'.format(len(file_data))
# Extract coordinates and zip them into lists.
x, y, mag = [], [], []
for st in file_data:
    # filter by mag.
    if st[k] < 20.:
        x.append(st[i])
        y.append(st[j])
        mag.append(st[k])
print 'Stars filtered by mag: {}'.format(len(x))
ra_cent = (max(x) + min(x)) / 2.
dec_cent = (max(y) + min(y)) / 2.
print 'Center of frame: RA={}, DEC={}'.format(ra_cent, dec_cent)
# Move coordinates to origin in center of frame. Apply cos(Dec) correction to
# RA coordinates.
x = (np.array(x) - ra_cent) * np.cos(np.deg2rad(dec_cent))
y = np.array(y) - dec_cent

# # Load Voronoi centers and radii.
# vor_file = 'out_fig_dat/junk_17.5_0.45_1.0_5.out'
# file_data = np.loadtxt(vor_file)
i, j, k = 4, 3, 5
# Extract coordinates and zip them into lists.
# c_x, c_y, r = zip(*file_data)[i], zip(*file_data)[j], zip(*file_data)[k]
c_x, c_y, r = [0.], [0.], [0.]

# Load BIca centers and radii.
com_file = '../cent_rads.dat'
file_data = np.loadtxt(com_file)
i, j, k, q = 0, 1, 2, 3
# Extract coordinates and zip them into lists.
cc_x, cc_y, cr_x, cr_y = zip(*file_data)[i], zip(*file_data)[j], \
    zip(*file_data)[k], zip(*file_data)[q]
# Move center coordinates to new origin.
cc_x = (np.array(cc_x) - ra_cent) * np.cos(np.deg2rad(dec_cent))
cc_y = (np.array(cc_y) - dec_cent)
# Convert RA, DEC longitudes in arcmin to radii in degrees.
cr_x, cr_y = np.array(cr_x) / (60. * 2), np.array(cr_y) / (60. * 2)

# Make plot.
data1, data2 = zip(*[c_x, c_y, r]), zip(*[cc_x, cc_y, cr_x, cr_y])
make_plot(x, y, mag, data1, data2)

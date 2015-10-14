

def save_cent_rad(f_name, a_f, cent_rad):
    '''
    Save center and radius data to file for each parameter value processed.
    '''

    data_out = 'out_fig_dat/' + f_name + '_' + str(a_f) + '.out'
    with open(data_out, 'w') as f:
        f.write("#area_frac   min_neigh    cent_x    cent_y    rad\n")
    with open(data_out, 'a') as f:
        for i, l in enumerate(cent_rad[0]):
            m_n = cent_rad[1][i]
            f.write("{:.3f}{:>11.0f}{:>16.4f}{:>10.4f}{:>10.4f}\n".format(
                a_f, m_n, *l))


def save_to_log(f_name, text, i):
    '''
    Save info to .log file.
    '''
    data_out = 'out_fig_dat/' + f_name + '.log'
    if i == 0:
        with open(data_out, 'w') as f:
            f.write("{}\n".format(text))
    else:
        with open(data_out, 'a') as f:
            f.write("{}\n".format(text))

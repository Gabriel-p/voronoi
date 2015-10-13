

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

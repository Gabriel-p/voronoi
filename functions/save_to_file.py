

def save_cent_rad(f_name, m_rang, area_frac_range, m_n, new_cent_rad,
                  old_cent_rad):
    '''
    Save center and radius data to file for each parameter value processed.
    '''

    a_f_min, a_f_max = [round(_, 2) for _ in area_frac_range]
    data_out = 'out_fig_dat/' + f_name + '_' + str(round(m_rang, 1)) + '_' +\
        str(a_f_min) + '_' + str(a_f_max) + '_' + str(m_n) + '.out'
    with open(data_out, 'w') as f:
        f.write("#a_f_min    a_f_max  min_n    cent_x    cent_y     rad\n")
    with open(data_out, 'a') as f:
        for i, l in enumerate(new_cent_rad):
            f.write("{:8.2f}{:8.2f}{:>6.0f}{:>16.4f}{:>10.4f}{:>10.4f}\n".
                    format(a_f_min, a_f_max, m_n, *l))

    if old_cent_rad:
        with open(data_out, 'a') as f:
            f.write("#REJECTED GROUPS BELOW.\n")
        with open(data_out, 'a') as f:
            for i, l in enumerate(old_cent_rad):
                f.write("{:8.2f}{:8.2f}{:>6.0f}{:>16.4f}{:>10.4f}{:>10.4f}\n".
                        format(a_f_min, a_f_max, m_n, *l))


def save_to_log(f_name, text, m_rang, i):
    '''
    Save info to .log file.
    '''
    data_out = 'out_fig_dat/' + f_name + '_' + str(round(m_rang, 1)) + '.log'
    if i == 0:
        with open(data_out, 'w') as f:
            f.write("{}\n".format(text))
    else:
        with open(data_out, 'a') as f:
            f.write("{}\n".format(text))

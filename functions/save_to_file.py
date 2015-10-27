

def save_cent_rad(f_name, cent_rad, dens_accp_groups, dens_rej_groups,
                  intens_accp_groups, intens_rej_groups):
    '''
    Save center and radius data to file for each parameter value processed.
    '''

    data_out = 'out_fig_dat/' + f_name + '.out'
    with open(data_out, 'w') as f:
        f.write("#\n#cent_x    cent_y     rad\n#\n")

    with open(data_out, 'a') as f:
        f.write("#All groups detected.\n")
        for l in cent_rad:
            f.write("{:<10.4f}{:>10.4f}{:>10.4f}\n".format(*l))

        f.write("#\n#Groups with accepted density values.\n")
        for l in dens_accp_groups:
            f.write("{:<10.4f}{:>10.4f}{:>10.4f}\n".format(*l))

        f.write("#\n#Groups with rejected density values.\n")
        for l in dens_rej_groups:
            f.write("{:<10.4f}{:>10.4f}{:>10.4f}\n".format(*l))

        f.write("#\n#Groups with accepted intensity/area values.\n")
        for l in intens_accp_groups:
            f.write("{:<10.4f}{:>10.4f}{:>10.4f}\n".format(*l))

        f.write("#\n#Groups with rejected intensity/area values.\n")
        for l in intens_rej_groups:
            f.write("{:<10.4f}{:>10.4f}{:>10.4f}\n".format(*l))


def save_to_log(f_name, text, w_a):
    '''
    Save info to .log file.
    '''
    print text
    data_out = 'out_fig_dat/' + f_name + '.log'
    with open(data_out, w_a) as f:
        f.write("{}\n".format(text))

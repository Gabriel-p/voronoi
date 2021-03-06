

def save_cent_rad(f_name, cent_rad, intens_acc_dens_acc, intens_acc_dens_rej,
                  intens_rej_dens_acc, intens_rej_dens_rej):
    '''
    Save center and radius data to file for each parameter value processed.
    '''

    data_out = 'out_fig_dat/' + f_name + '.out'
    with open(data_out, 'w') as f:
        f.write("#\n")

    with open(data_out, 'a') as f:
        f.write("# All groups detected.\n")
        f.write("#\n#     cent_x         cent_y            rad\n#\n")
        for l in cent_rad:
            f.write("{:>12.6f}{:>15.6f}{:>15.6f}\n".format(*l))

        id_line = "#\n#     cent_x         cent_y            rad       N/A" +\
            "       I/A\n#\n"
        l_format = "{:>12.6f}{:>15.6f}{:>15.6f}{:>10.4f}{:>10.4f}\n"

        f.write("#\n# Groups with accepted density + I/A values.\n")
        f.write(id_line)
        for i, l in enumerate(intens_acc_dens_acc[0]):
            l += [intens_acc_dens_acc[2][i], intens_acc_dens_acc[3][i]]
            f.write(l_format.format(*l))

        f.write("#\n# Groups with accepted density + rejected I/A values.\n")
        f.write(id_line)
        for i, l in enumerate(intens_rej_dens_acc[0]):
            l += [intens_rej_dens_acc[2][i], intens_rej_dens_acc[3][i]]
            f.write(l_format.format(*l))

        f.write("#\n# Groups with rejected density + accepted I/A values.\n")
        f.write(id_line)
        for i, l in enumerate(intens_acc_dens_rej[0]):
            l += [intens_acc_dens_rej[2][i], intens_acc_dens_rej[3][i]]
            f.write(l_format.format(*l))

        f.write("#\n# Groups with rejected density + I/A values.\n")
        f.write(id_line)
        for i, l in enumerate(intens_rej_dens_rej[0]):
            l += [intens_rej_dens_rej[2][i], intens_rej_dens_rej[3][i]]
            f.write(l_format.format(*l))


def save_to_log(f_name, text, w_a):
    '''
    Save info to .log file.
    '''
    print text
    data_out = 'out_fig_dat/' + f_name + '.log'
    with open(data_out, w_a) as f:
        f.write("{}\n".format(text))

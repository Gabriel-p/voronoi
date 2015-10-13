
import re
import numpy as np


def char_remove(in_lst):
    '''
    Correctly convert input data for parameters ranges to lists.
    '''
    lst = []
    # If input list is empty, return empty list.
    if in_lst[1:]:
        l0 = []
        if in_lst[1][0] in {'[', '(', '{'}:
            # Remove non-numeric characters and append numbers as floats.
            l0.append([float(i) for i in re.findall('[0-9.]+',
                      str(in_lst[1:]))])
            # Store indicating that this is a list of values.
            lst = ['l', map(float, l0[0])]
        else:
            # Store indicating that this is a range of values.
            lst = ['r', map(float, in_lst[1:4])]

    return lst


def get_vals(params):
    '''
    Convert lists to proper values if their are ranges, or return clean
    lists if they are already lists.
    '''
    # If it's a range, obtain that range.
    if params[0] == 'r':
        par_vals = np.arange(params[1][0], params[1][1], params[1][2])
    else:
        par_vals = params[1]

    return par_vals


def get_params_in():

    # Read data from file.
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):

            if not line.startswith("#") and line.strip() != '':
                reader = line.split()

                # Updater.
                if reader[0] == 'MA':
                    avr_area_frac = char_remove(reader)

                # Mode.
                elif reader[0] == 'MN':
                    min_neighbors = char_remove(reader)

    avr_area_frac = get_vals(avr_area_frac)
    min_neighbors = get_vals(min_neighbors)

    return avr_area_frac, min_neighbors

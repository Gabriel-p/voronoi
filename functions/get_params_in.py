

def get_params_in():

    # Read data from file.
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):

            if not line.startswith("#") and line.strip() != '':
                reader = line.split()

                # Maximum area.
                if reader[0] == 'FN':
                    in_file = str(reader[1])

                # Magnitude range.
                elif reader[0] == 'MR':
                    mag_range = map(float, reader[1:])

                # Area range.
                elif reader[0] == 'MA':
                    area_frac_range = map(float, reader[1:])

                # Minimum neighbors.
                elif reader[0] == 'MN':
                    min_neighbors = int(reader[1])

                # Minimum neighbors.
                elif reader[0] == 'FI':
                    intens_frac = float(reader[1])

    return in_file, mag_range, area_frac_range, min_neighbors, intens_frac

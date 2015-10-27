

def get_params_in():

    # Read data from file.
    with open('params_input.dat', "r") as f_dat:
        # Iterate through each line in the file.
        for l, line in enumerate(f_dat):

            if not line.startswith("#") and line.strip() != '':
                reader = line.split()

                # Input file's name.
                if reader[0] == 'FN':
                    in_file = str(reader[1])
                # Input file's data columns.
                if reader[0] == 'DC':
                    in_file_cols = map(int, reader[1:-1])
                    coords_flag = str(reader[-1])

                # Magnitude range.
                elif reader[0] == 'MR':
                    mag_range = map(float, reader[1:])

                # Area range.
                elif reader[0] == 'MA':
                    area_frac_range = map(float, reader[1:])

                # Minimum neighbors.
                elif reader[0] == 'MN':
                    min_neighbors = int(reader[1])

                # Fraction of frame's intensity/area.
                elif reader[0] == 'FI':
                    intens_frac = float(reader[1])

                # Fraction of frame's density (stars/area).
                elif reader[0] == 'FD':
                    dens_frac = float(reader[1])

    return in_file, in_file_cols, coords_flag, mag_range, area_frac_range,\
        min_neighbors, intens_frac, dens_frac

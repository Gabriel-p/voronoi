

def print_perc(i, tot_sols, milestones):
    '''
    '''
    percentage_complete = (100.0 * (i + 1) / tot_sols)
    while len(milestones) > 0 and \
            percentage_complete >= milestones[0]:
        print "    {:>3}%".format(milestones[0])
        # Remove that milestone from the list.
        milestones = milestones[1:]

    return milestones

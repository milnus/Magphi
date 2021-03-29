import os


def check_number_of_primers(primer_file):
    # Initiate primer counter
    num_primers = 0

    # Go through primer file
    with open(primer_file, 'r') as primers:
        for line in primers.readlines():
            # Check if a new primer has been reached
            if ">" in line:
                num_primers += 1

    if num_primers % 2 != 0:
        print("UPS! The number of primers given is not even! This means that not all primers can be given a mate\n"
              "If you have a primer that should be used twice it should be in the file twice.")
        exit(code=1)
    else:
        return num_primers


def extract_primer_info(primer_file):
    # Initiate dict to hold primers with their name as the key and sequence as the value.
    primer_dict = {}

    # Go through primer file
    with open(primer_file, 'r') as primers:
        for line in primers.readlines():
            # Remove new lines from line
            line = line.strip('\n')
            # Check if a new primer has been reached, if not then add line of sequence to current primer
            if ">" in line:
                primer_name = line[1:]
                primer_dict[primer_name] = []
            else:
                # TODO - AT THE MOMENT WE DO NOT USE THE SEQUENCE OF THE PRIMER AS WE CAN BLAST ALL PRIMERS AT ONCE. - Maybe keep them if we want to do the iterative decrease of blast parameters to get a specific primer pair to fit.
                primer_dict[primer_name].append(line)

    return primer_dict


def construct_pair_primers(primer_names):
    # Initiate dict that contain the basename of a primer pair as the key and primer names as the values
    primer_pairs = {}

    while len(primer_names) > 0:
        basename_lengths = []
        chosen_primer = primer_names.pop(0)

        # Go through primers and find the length of the base name compared to the chosen
        for i, primer in enumerate(primer_names):
            # compare chosen primer to a possible mate
            basename_lengths.append(len(os.path.commonprefix([chosen_primer, primer])))

        # Find the primer(s) with the longest base name
        longest_mates_index = [mate_index for mate_index, mate in enumerate(basename_lengths) if mate == max(basename_lengths)]

        # Check if only one mate is found, if then record the pair, if not then return the chosen one to the list a restart
        if len(longest_mates_index) == 1:
            # Get the best match
            chosen_mate = primer_names[longest_mates_index[0]]

            # Construct the pair
            primer_pair = [chosen_mate, chosen_primer]

            # Find the common name for pair
            common_name = os.path.commonprefix(primer_pair)

            if common_name[-1] == '_':
                common_name = common_name[:-1]

            # Insert pair into primer dict
            primer_pairs[common_name] = primer_pair

            # Delete primer chosen as mate from list of primers
            del primer_names[longest_mates_index[0]]
        else:
            primer_names.append(chosen_primer)

    return primer_pairs


def handle_primers(primer_file):
    num_primers = check_number_of_primers(primer_file)

    primer_dict = extract_primer_info(primer_file)

    primer_pairs = construct_pair_primers(list(primer_dict.keys()))

    # TODO - make verbose controlled:
    print(f'{num_primers} primers found in primer file.\n'
          f'This leads to {len(primer_pairs)} pairs of primers.')

    return primer_pairs, primer_dict

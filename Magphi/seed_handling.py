import os
import sys
try:
    from Magphi.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error

EXIT_INPUT_FILE_ERROR = 3


def check_number_n_names_of_seeds(seed_file, file_logger):
    """
    Function to count and report on the number of seed sequences and pairs of seeds expected
    Also record all names of seeds sequences in a list
    :param
    seed_file: File name of multi fasta file containing seed sequences
    file_logger: Logger that outputs files to log
    :return: Will throw error is number of seed sequences is odd
    or if two seeds have an identical name,
    else returns the number of seeds and the names in a list
    """
    file_logger.debug("Initiating count of seed sequences")
    # Initiate seed counter
    num_seeds = 0
    seed_list = []

    # Go through multi fasta file of seed sequences
    with open(seed_file, 'r') as seeds:
        for line in seeds.readlines():
            # Check if a new seed sequence has been reached, if then increment seed sequence count
            # and add the name to the list of seeds
            if ">" in line:
                num_seeds += 1

                line = line.strip('\n')
                seed_name = line[1:]

                if seed_name in seed_list:
                    file_logger.exception(
                        "Duplicate seed sequence names were identified! This is not allowed and should be changed")
                    exit_with_error(
                        message="Duplicate seed sequence names were identified! This is not allowed and should be changed",
                        exit_status=EXIT_INPUT_FILE_ERROR)

                seed_list.append(seed_name)

    if num_seeds % 2 != 0:
        file_logger.exception("The number of seed sequences given is not even! This means that not all seed sequences can be given a mate\n"
              "If you have a seed sequence that should be used twice it should be in the file twice.")
        exit_with_error(message="The number of seed sequences given is not even! This means that not all seed sequences can be given a mate\n"
                                "If you have a seed sequence that should be used twice it should be in the file twice.",
                                        exit_status=EXIT_INPUT_FILE_ERROR)
    else:
        file_logger.debug("Number of seed sequences given is equal")
        return num_seeds, seed_list


def construct_pair_seeds(seed_names, file_logger):
    """
    Function to pair names of seed sequences into pairs of two.
    :param seed_names: Names of seed sequences from input file
    :param file_logger: Logger that outputs files to log
    :return: A dict with key names being seed sequence pair name and values the name of seed sequences
    """
    file_logger.debug("Initiate pairing of seed sequence names into pairs")
    # Initiate dict that contain the basename of a seed pair as the key and seed names as the values
    seed_pairs = {}

    # Loop through all seed sequence names to find pairs of seed sequences
    n_loops = 0
    while len(seed_names) > 0 and n_loops != 1000:
        basename_lengths = []
        # Choose and extract fist seed in list
        chosen_seed = seed_names.pop(0)

        # Go through seed sequences and find the length of the base name compared to the chosen
        for i, seed in enumerate(seed_names):
            # Compare chosen seed to a possible mate and extract the length of base/common name
            basename_lengths.append(len(os.path.commonprefix([chosen_seed, seed])))

        # Find the seed(s) with the longest base name
        longest_mates_index = [mate_index for mate_index, mate in enumerate(basename_lengths) if mate == max(basename_lengths)]

        # Check if only one seed sequence is found to have the longest base/common name with chosen seed sequence,
        #   if then record the pair,
        #   if not then return the chosen one to the list a restart
        if len(longest_mates_index) == 1:
            # Get the best match
            chosen_mate = seed_names[longest_mates_index[0]]

            # Construct the pair
            seed_pair = [chosen_mate, chosen_seed]

            # Find the common name for pair
            common_name = os.path.commonprefix(seed_pair)

            if common_name[-1] == '_':
                common_name = common_name[:-1]

            # Insert pair into seed dict
            seed_pairs[common_name] = seed_pair

            # Delete seed chosen as mate from list of seeds
            del seed_names[longest_mates_index[0]]
        else:
            seed_names.append(chosen_seed)

        n_loops += 1

    if n_loops == 1000:
        file_logger.exception(f'Seed sequence pairing failed due to more than 1000 rounds of pairing tested. seeds remaining to be paired: seed_names = {seed_names}')
        exit_with_error(message=f'Seed sequence pairing failed due to more than 1000 rounds of pairing tested. seeds remaining to be paired: seed_names = {seed_names}',
                                        exit_status=EXIT_INPUT_FILE_ERROR)

    return seed_pairs


def handle_seeds(seed_file, file_logger):
    """
    Function to be called in main that will check, extract info from, and pair seed sequences into pairs
    returns the Seed sequence pairs and a dict with name and sequences of seeds
    :param seed_file: File name of multi sequence fasta file containing the seed sequences
    :param file_logger: Logger that outputs files to log
    :return: pairs of seed sequences
    """

    num_seeds, seed_list = check_number_n_names_of_seeds(seed_file, file_logger)

    seed_pairs = construct_pair_seeds(seed_list, file_logger)

    file_logger.debug(f'{num_seeds} seeds found in seed file. - '
                      f'This leads to {len(seed_pairs)} pairs of seeds.')

    return seed_pairs

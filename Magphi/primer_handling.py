import os
import sys
try:
    from Magphi.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error

EXIT_INPUT_FILE_ERROR = 3


def check_number_n_names_of_seeds(primer_file, file_logger):
    """
    Function to count and report on the number of seed sequences and pairs of primers expected
    Also record all names of seeds sequences in a list
    :param
    primer_file: File name of multi fasta file containing seed sequences
    file_logger: Logger that outputs files to log
    :return: Will throw error is number of seed sequences is odd
    or if two seeds have an identical name,
    else returns the number of primers and the names in a list
    """
    file_logger.debug("Initiating count of seed sequences")
    # Initiate primer counter
    num_primers = 0
    seed_list = []

    # Go through multi fasta file of seed sequences
    with open(primer_file, 'r') as primers:
        for line in primers.readlines():
            # Check if a new seed sequence has been reached, if then increment seed sequence count
            # and add the name to the list of seeds
            if ">" in line:
                num_primers += 1

                line = line.strip('\n')
                primer_name = line[1:]

                if primer_name in seed_list:
                    file_logger.exception(
                        "Duplicate seed sequence names were identified! This is not allowed and should be changed")
                    exit_with_error(
                        message="Duplicate seed sequence names were identified! This is not allowed and should be changed",
                        exit_status=EXIT_INPUT_FILE_ERROR)

                seed_list.append(primer_name)

    if num_primers % 2 != 0:
        file_logger.exception("The number of seed sequences given is not even! This means that not all seed sequences can be given a mate\n"
              "If you have a seed sequence that should be used twice it should be in the file twice.")
        exit_with_error(message="The number of seed sequences given is not even! This means that not all seed sequences can be given a mate\n"
                                "If you have a seed sequence that should be used twice it should be in the file twice.",
                                        exit_status=EXIT_INPUT_FILE_ERROR)
    else:
        file_logger.debug("Number of seed sequences given is equal")
        return num_primers, seed_list


def construct_pair_primers(primer_names, file_logger):
    """
    Function to pair names of seed sequences into pairs of two.
    :param primer_names: Names of seed sequences from input file
    :param file_logger: Logger that outputs files to log
    :return: A dict with key names being seed sequence pair name and values the name of seed sequences
    """
    file_logger.debug("Initiate pairing of seed sequence names into pairs")
    # Initiate dict that contain the basename of a primer pair as the key and primer names as the values
    primer_pairs = {}

    # Loop through all seed sequence names to find pairs of seed sequences
    n_loops = 0
    while len(primer_names) > 0 and n_loops != 1000:
        basename_lengths = []
        # Choose and extract fist primer in list
        chosen_primer = primer_names.pop(0)

        # Go through seed sequences and find the length of the base name compared to the chosen
        for i, primer in enumerate(primer_names):
            # Compare chosen primer to a possible mate and extract the length of base/common name
            basename_lengths.append(len(os.path.commonprefix([chosen_primer, primer])))

        # Find the primer(s) with the longest base name
        longest_mates_index = [mate_index for mate_index, mate in enumerate(basename_lengths) if mate == max(basename_lengths)]

        # Check if only one seed sequence is found to have the longest base/common name with chosen seed sequence,
        #   if then record the pair,
        #   if not then return the chosen one to the list a restart
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

        n_loops += 1

    if n_loops == 1000:
        file_logger.exception(f'Seed sequence pairing failed due to more than 1000 rounds of pairing tested. primers remaining to be paired: primer_names = {primer_names}')
        exit_with_error(message=f'Seed sequence pairing failed due to more than 1000 rounds of pairing tested. primers remaining to be paired: primer_names = {primer_names}',
                                        exit_status=EXIT_INPUT_FILE_ERROR)

    return primer_pairs


def handle_primers(primer_file, file_logger):
    """
    Function to be called in main that will check, extract info from, and pair seed sequences into pairs
    returns the Seed sequence pairs and a dict with name and sequences of primers
    :param primer_file: File name of multi sequence fasta file containing the seed sequences
    :param file_logger: Logger that outputs files to log
    :return: pairs of seed sequences
    """

    num_primers, seed_list = check_number_n_names_of_seeds(primer_file, file_logger)

    primer_pairs = construct_pair_primers(seed_list, file_logger)

    file_logger.debug(f'{num_primers} primers found in primer file. - '
                      f'This leads to {len(primer_pairs)} pairs of primers.')

    return primer_pairs

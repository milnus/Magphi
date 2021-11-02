import gzip
try:
    from Magphi.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error

EXIT_INPUT_FILE_ERROR = 3


def check_if_gzip(input_files, file_logger):
    """
    Check if genome files given as input at gzipped, by trying to open them with zgip's open() funciton
    :param
    input_files: A list of names of genomes given as input.
    file_logger: A logger that will write to a logging file
    :return: Boolean that tell if all input files are gzipped.
    """
    file_logger.debug('Check if input genomes are in gzipped state')
    # Make list with Bool for each input genome and evaluate all input genomes for Gzipped state,
    # If genome genome is Gzipped change list index for genome to true in is_input_gzipepd.
    is_input_gzipped = [False] * len(input_files)
    for i, file in enumerate(input_files):
        with gzip.open(file, 'r') as input_file:
            try:
                input_file.readline()
                is_input_gzipped[i] = True
            except OSError:
            # except gzip.BadGzipFile:
                is_input_gzipped[i] = False

    # Check if all inputs are gzipped,
    # if then return True for all inputs being gzipped.
    # Else check if some are gzipped and return that input is mixed
    # Lastly return False if no input is Gzipped.
    if all(is_input_gzipped):
        file_logger.debug('Files were found to be gzipped')
        return True
    elif 0 < sum(is_input_gzipped) < len(is_input_gzipped):
        file_logger.exception('Check input files: Some input files found to be gzipped!\n'
                          'Please check that all of your input files are either Gzipped or not')
        exit_with_error(message='Some input files found to be Gzipped while other were not!',
                        exit_status=EXIT_INPUT_FILE_ERROR)
    else:
        file_logger.debug('Files were not found to be gzipped')
        return False


def check_if_fasta(input_files, file_logger, is_input_gzipped):
    """
    Function to check if input files are Fasta format identified by > in first line
    returns 'fasta' if all files are fasta, exits and logs error if only some are fasta indicating mixed input
    :param input_files: List of input genomes by file name
    :param file_logger: Logger that outputs files to log
    :return: variable: Either 'fasta' if all input genomes are fasta format or None if not
    """
    file_logger.debug('Check if input genomes are fasta format')

    # Construct list to hold if files are fasta
    # Search all files for fasta signature > in first line.
    is_input_fasta = [False]*len(input_files)
    for i, file in enumerate(input_files):
        if is_input_gzipped:
            in_file = gzip.open(file, 'rt')
        else:
            in_file = open(file, 'r')

        # Test first line
        first_line = in_file.readline()
        if '>' in first_line:
            if len(first_line.strip()) > 1:
                is_input_fasta[i] = True
            else:
                in_file.close()
                exit_with_error(
                    f'Fasta file contains no sequence header. This is not allowed please have a look at file: {file}',
                    EXIT_INPUT_FILE_ERROR)

        # Check for new line in remaining lines
        new_line_in_middle = False
        for line in in_file.readlines():
            if new_line_in_middle:
                in_file.close()
                exit_with_error(f'Fasta file contains new line in middle of file. This is not allowed please have a look at file: {file}',
                                EXIT_INPUT_FILE_ERROR)
            # Check for empty line
            if not line.strip():
                new_line_in_middle = True
        in_file.close()


    # Check if all input genomes are fasta, if then return 'fasta'
    # else check if only some input genomes are fasta - if: Give error
    # Else return None
    if all(is_input_fasta):
        file_logger.debug('Files were found to be Fasta files')
        return 'fasta'

    elif any(is_input_fasta):
        file_logger.error('Input files were found to be of mixed type with some being recognised as fasta!\n'
                      'Please make sure that all your input files are only either Fasta or GFF3, not mixed')
        exit_with_error(message='Input files found to be of mixed type! Only some files were recognised as fasta!',
                        exit_status=EXIT_INPUT_FILE_ERROR)

    else:
        file_logger.debug('Files were not found to be Fasta files')
        return None


def check_if_gff(input_files, file_logger, is_input_gzipped):
    """
    Function to check in input genomes are given in a gff3 format with the appended genome.
    :param
    input_files: List of input genomes by file name
    file_logger: Logger that outputs files to log
    :return: Either 'gf' if all input genomes are fasta format or None if not
    """
    ''' Function to check if input files are GFF3 format with an attached genomes identified by the ##FASTA line '''
    file_logger.debug('Check input files: Check if GFF3')
    is_input_gff = [False] * len(input_files)
    for i, file in enumerate(input_files):
        if is_input_gzipped:
            in_file = gzip.open(file, 'rt')
        else:
            in_file = open(file, 'r')

        if '##gff-version 3' in in_file.readline():
            for line in in_file.readlines():
                if '##FASTA' in line:
                    is_input_gff[i] = True
            # Check that a genome has been found in file
            if is_input_gff[i] is False:
                in_file.close()
                print(f'UPS! {file} does seem to be a GFF3 version, but not one that contain the genome '
                      f'following a ##FASTA line - Please check the file')
                file_logger.error('Check input files: Input is GFF but seems to miss ##FASTA line or entire genome')
                exit_with_error(message='Input is GFF but seems to miss ##FASTA line or entire genome',
                                exit_status=EXIT_INPUT_FILE_ERROR)
        in_file.close()

    if all(is_input_gff):
        file_logger.debug('Files were found to be GFF files')
        return 'gff'
    elif any(is_input_gff):
        file_logger.error('Input files were found to be of mixed type with some being recognised as GFF!\n'
                      'Please make sure that all your input files are only either Fasta or GFF3, not mixed')
        exit_with_error(message='Input files found to be of mixed type! Only some files were recognised as GFF!',
                        exit_status=EXIT_INPUT_FILE_ERROR)
    else:
        file_logger.debug('Files were not found to be GFF files')
        return None


def check_inputs(input_files, file_logger):
    """
    Function to run through different checks of input files to determine their state (gzipped) and type (Fasta or GFF)
    :param input_files: List of input genomes by file name
    :param file_logger: Logger that outputs files to log
    :return: Variable for Type, and state of file found by checks.
    """
    ''' Function to check the input files. Will chack if input is GFF, Fasta or exit if filetype is not recognised.
     Will exit program when unrecognised input file is given.'''

    # Add logging step
    file_logger.debug("Check input files: "
                      "Initiating verifying input files are either FASTA or GFF3 with genome attached."
                      " Also checking for gzipped state")

    # Check if files are gzipped
    is_input_gzipped = check_if_gzip(input_files, file_logger)

    # check if input file is fasta
    file_type = check_if_fasta(input_files, file_logger, is_input_gzipped)

    # Check if input files are GFF3
    if file_type is None:
        file_type = check_if_gff(input_files, file_logger, is_input_gzipped)

    # If file_type is still None, exit as input files are not recognized
    if file_type is None:
        file_logger.exception('The given input files could not be recognised as either Fasta or GFF3 with a genome.')
        exit_with_error(message='The given input files could not be recognised as either Fasta or GFF3 with a genome',
                        exit_status=EXIT_INPUT_FILE_ERROR)

    return file_type, is_input_gzipped

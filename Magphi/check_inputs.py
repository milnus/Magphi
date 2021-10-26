import logging
import gzip
try:
    from Magphi.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error

EXIT_INPUT_FILE_ERROR = 3


def check_if_gzip(input_files):
    is_input_gzipped = [False] * len(input_files)
    for i, file in enumerate(input_files):
        with gzip.open(file, 'r') as input_file:
            try:
                input_file.readline()
                is_input_gzipped[i] = True
            except OSError:
            # except gzip.BadGzipFile:
                is_input_gzipped[i] = False

    if all(is_input_gzipped):
        return True
    elif 0 < sum(is_input_gzipped) < len(is_input_gzipped):
        print('UPS! Something went wrong!\n'
              'Some input files found to be gzipped! '
              'Please check that all of your input files are either Gzipped or not')
        logging.error('Check input files: Some input files found to be gzipped!\n'
                      'Please check that all of your input files are either Gzipped or not')
        exit_with_error(message='Some input files found to be Gzipped while other were not!',
                        exit_status=EXIT_INPUT_FILE_ERROR)
    else:
        return False


def check_if_fasta(input_files):
    ''' Function to check if input files are Fasta format identified by > in first line
    returns 'fasta' if all files are fasta, exits and logs error if only some are fasta indicating mixed input '''
    logging.info('Check input files: Check if fasta')
    # Construct list to hold if files are fasta
    is_input_fasta = [False]*len(input_files)

    # Go through each file and check if it is a fasta genome
    for i, file in enumerate(input_files):
        with open(file, 'r') as in_file:
            if '>' in in_file.readline():
                is_input_fasta[i] = True

    if all(is_input_fasta):
        # TODO - make verbose controlled
        print('Input files found to be Fasta format')
        logging.info('Check input files: Files are found to be Fasta files')
        return 'fasta'

    elif any(is_input_fasta):
        print('UPS! Something went wrong!\n'
              'Input files found to be of mixed type! '
              'Please check that all your input files are either Fasta or GFF, not mixed')
        logging.error('Check input files: Input files found to be of mixed type!\n'
                      'Please make sure that all your input files are either Fasta or GFF3, not mixed')
        exit_with_error(message='Input files found to be of mixed type!',
                        exit_status=EXIT_INPUT_FILE_ERROR)

    else:
        return None


def check_if_gff(input_files):
    ''' Function to check if input files are GFF3 format with an attached genomes identified by the ##FASTA line '''
    logging.info('Check input files: Check if GFF3')
    is_input_gff = [False] * len(input_files)
    for i, file in enumerate(input_files):
        with open(file, 'r') as in_file:
            if '##gff-version 3' in in_file.readline():
                for line in in_file.readlines():
                    if '##FASTA' in line:
                        is_input_gff[i] = True
                # Check that a genome has been found in file
                if is_input_gff[i] is False:
                    print(f'UPS! {file} does seem to be a GFF3 version, but not one that contain the genome '
                          f'following a ##FASTA line - Please check the file')
                    logging.error('Check input files: Input is GFF but seems to miss ##FASTA line or entire genome')
                    exit_with_error(message='Input is GFF but seems to miss ##FASTA line or entire genome',
                                    exit_status=EXIT_INPUT_FILE_ERROR)

    if all(is_input_gff):
        # TODO - make verbose controlled
        print(f'all files were found to be GFF3 files')
        logging.info('Check input files: Files are found to be GFF3 files')
        return 'gff'
    elif any(is_input_gff):
        print('UPS! Something went wrong!\n'
              'Some files seems to not satisfy GFF3 with attached genome! '
              'Please check that all your input files are either Fasta or GFF, not mixed')
        logging.error('Check input files: Input files found to be of mixed type!\n'
                      'Please make sure that all your input files are either Fasta or GFF3, not mixed')
        exit_with_error(message='Input files found to be of mixed type some GFF some not!',
                        exit_status=EXIT_INPUT_FILE_ERROR)
    else:
        return None


def check_inputs(input_files):
    ''' Function to check the input files. Will chack if input is GFF, Fasta or exit if filetype is not recognised.
     Will exit program when unrecognised input file is given.'''

    # Add logging step
    logging.info("Check input files: Initiating verifying input files are either FASTA or GFF3 with genome attached")

    # Check if files are gzipped
    is_input_gzipped = check_if_gzip(input_files)

    # check if input file is fasta
    file_type = check_if_fasta(input_files)

    # Check if input files are GFF3
    if file_type is None:
        file_type = check_if_gff(input_files)

    # If file_type is still none, exit as input files are not recognized
    if file_type is None:
        logging.error('Input check: The given input files could not be recognised as either Fasta or GFF3 with a genome.')
        exit_with_error(message='The given input files could not be recognised as either Fasta or GFF3 with a genome',
                        exit_status=EXIT_INPUT_FILE_ERROR)

    return file_type, is_input_gzipped

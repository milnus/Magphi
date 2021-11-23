import argparse
import sys
try:
    from Magphi.check_depencies import check_dependencies_only
except ModuleNotFoundError:
    from check_depencies import check_dependencies_only

EXIT_COMMAND_LINE_ERROR = 2

# def get_commandline_arguments(args, version):
def get_commandline_arguments(args, version):
    ''' Parse command line arguments.
    return instructions given on command line in the args argument.
    Will recognise the unusual '-help' can and print help if nothing is given as input on the command line.'''
    # Set up parser
    parser = argparse.ArgumentParser(prog='Magphi',
                                     description='Welcome to Magphi!\n '
                                                 'This program will extract sequences and possible annotations within '
                                                 'a set of seed sequences given as input.')

    # Add the flag for the input genomes
    parser.add_argument('-g',
                        '--input_genomes',
                        help='Give the fasta or gff3 files. (gff3 files should contain the genome fasta sequence)',
                        dest='genomes',
                        required=True,
                        nargs='+',
                        metavar='.fa/.gff')

    # Add the flag for the multi fasta file containing seed sequences
    parser.add_argument('-s',
                        '--input_seeds',
                        help='Give the multi fasta containing the seed sequences to be used for extracting sequnces',
                        dest='seeds',
                        required=True,
                        metavar='multi_fasta_file.fa')

    parser.add_argument('-is',
                        '--include_seeds',
                        help='Argument to include the seeds in the sequence/annotations extracted '
                             '[default: seeds are removed]',
                        dest='include_seeds',
                        required=False,
                        action='store_true',
                        default=False)

    parser.add_argument('-md',
                        '--max_seed_distance',
                        help='The maximum distance with which seeds will be merged.\n'
                             'This can often be set a bit higher than an expected size of a region\n'
                             'If no maximum distance is wanted then set to 0 [default: 50,000bp]',
                        default=50000,
                        required=False,
                        type=int,
                        metavar='int',
                        dest='max_seed_dist')

    output_amount = parser.add_mutually_exclusive_group()
    output_amount.add_argument('-b',
                               '--print_breaks',
                               help='Argument to print outputs when seeds are next to contig breaks'
                                    ' [default: sequences are not printed]',
                               dest='print_breaks',
                               required=False,
                               action='store_true',
                               default=False)

    output_amount.add_argument('-n',
                               '--no_sequences',
                               help='Argument to not print outputs related to annotations or sequences found between seeds'
                                    ' [default: sequences are printed]',
                               dest='no_seqs',
                               required=False,
                               action='store_false',
                               default=True)

    # Add the flag for the output folder
    parser.add_argument('-o',
                        '--output_folder',
                        help='Give path to output folder [default: current folder]',
                        required=False,
                        metavar='path/to/output',
                        default='.',
                        dest='out_path')

    # Add the flag to control max CPUs
    parser.add_argument('-c',
                        '--cpu',
                        help='Give max number of CPUs [default: 1]',
                        required=False,
                        metavar='int',
                        default=1,
                        type=int,
                        dest='cpu')

    logger_level = parser.add_mutually_exclusive_group()
    logger_level.add_argument('-l',
                              '--log',
                              help='Record program progress in for debugging purpose',
                              action='store_true',
                              default=False,
                              required=False)

    logger_level.add_argument('-q',
                              '--quiet',
                              help='Only print warnings',
                              action='store_true',
                              default=False,
                              required=False)

    parser.add_argument('--check',
                              help='Check dependencies for Magphi and exit',
                              dest='dependency_check',
                              action='store_true',
                              default=False,
                              required=False)

    parser.add_argument('-v',
                        '--version',
                        action='version',
                        version=f'Magphi {version}')

    # Check if there are no arguments given or the user ask for the help message
    if len(args) < 1:
        parser.print_help()
        sys.exit(EXIT_COMMAND_LINE_ERROR)
    elif '-help' in args:
        parser.print_help()
        sys.exit(0)
    if '--check' in args:
        check_dependencies_only()

    args = parser.parse_args(args)

    return args

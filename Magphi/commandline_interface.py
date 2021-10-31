import argparse
import sys

EXIT_COMMAND_LINE_ERROR = 2

# TODO - Change all primer mentions to seed seqeunce
# def get_commandline_arguments(args, version):
def get_commandline_arguments(args, version):
    ''' Parse command line arguments.
    return instructions given on command line in the args argument.
    Will recognise the unusual '-help' can and print help if nothing is given as input on the command line.'''
    # Set up parser
    parser = argparse.ArgumentParser(prog='Magphi',
                                     description='Welcome to Magphi!\n '
                                                 'This program will extract sequences and possible annotations within '
                                                 'a set of insertion sequences or primers given as input.')

    # Add the flag for the input genomes
    parser.add_argument('-g',
                        '--input_genomes',
                        help='Give the fasta or gff3 files. (gff3 files should contain the genome fasta sequence)',
                        dest='genomes',
                        required=True,
                        nargs='+',
                        metavar='.fa/.gff')

    # Add the flag for the multi fasta file containing insertion sequences or primers.
    parser.add_argument('-s',
                        '--input_seeds',
                        help='Give the multi fasta containing the primers to be used for extracting sequnces',
                        dest='seeds',
                        required=True,
                        metavar='multi_fasta_file.fa')

    parser.add_argument('-ip',
                        '--include_primers',
                        help='Argument to include the primers in the sequence/annotations extracted '
                             '[default: primers are removed]',
                        dest='include_primers',
                        required=False,
                        action='store_true',
                        default=False)

    parser.add_argument('-md',
                        '--max_primer_distance',
                        help='The maximum distance with which primers will be merged.\n'
                             'This can often be set a bit higher than an expected size of a region\n'
                             'If no maximum distance is wanted then set to 0 [default: 50,000bp]',
                        default=50000,
                        required=False,
                        type=int,
                        metavar='int',
                        dest='max_primer_dist')

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

    parser.add_argument('-l',
                        '--log',
                        help='Record program progress in for debugging purpose',
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

    args = parser.parse_args(args)

    return args

import argparse


def get_commandline_arguments(args):
    # Set up parser
    parser = argparse.ArgumentParser(description='Welcome to Phupha!\n '
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
    parser.add_argument('-p',
                        '--input_primers',
                        help='Give the multi fasta containing the primers to be used for extracting sequnces',
                        dest='primers',
                        required=True,
                        metavar='multi_fasta_file')

    parser.add_argument('-ip',
                        '--include_primers',
                        help='Argument to include the primers in the sequence/annotations extracted '
                             '[default: primers are removed]',
                        dest='include_primers',
                        required=False,
                        action='store_false',
                        default=True)

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
    # TODO - See if the nargs works if you set the type expected to something (str)

    # Check if there are no arguments given or the user ask for the help message
    if len(args) < 1 or '-help' in args:
        parser.print_help()
        exit(code=1)

    args = parser.parse_args(args)

    return args

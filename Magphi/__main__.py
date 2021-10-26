'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Magnus Ganer Jespersen, 11 Oct 2021 
License     : MIT 
Maintainer  : magnus.ganer.j@gmail.com 
Portability : POSIX

The program reads one or more input FASTA files. For each file it computes a
variety of statistics, and then prints a summary of the statistics as output.
'''

import warnings
import os
import time
import logging
from sys import argv
import pkg_resources # ??
import concurrent.futures

try:
    from Magphi.commandline_interface import get_commandline_arguments
except ModuleNotFoundError:
    from commandline_interface import get_commandline_arguments

try:
    from Magphi.check_depencies import check_dependencies_for_main
except ModuleNotFoundError:
    from check_depencies import check_dependencies_for_main

try:
    from Magphi.exit_with_error import exit_with_error
except ModuleNotFoundError:
    from exit_with_error import exit_with_error

try:
    from Magphi.check_inputs import check_inputs
except ModuleNotFoundError:
    from check_inputs import check_inputs

try:
    from Magphi.split_gff_file import split_gff_files
except ModuleNotFoundError:
    from split_gff_file import split_gff_files

try:
    from Magphi.primer_handling import handle_primers
except ModuleNotFoundError:
    from primer_handling import handle_primers

try:
    from Magphi.search_insertion_sites import screen_genome_for_primers
except ModuleNotFoundError:
    from search_insertion_sites import screen_genome_for_primers

try:
    from Magphi.wrangle_outputs import partition_outputs, write_paired_primers
except ModuleNotFoundError:
    from wrangle_outputs import partition_outputs, write_paired_primers

try:
    from Magphi.write_output_csv import write_primer_hit_matrix, write_annotation_num_matrix, write_primer_hit_evidence, write_inter_primer_dist
except ModuleNotFoundError:
    from write_output_csv import write_primer_hit_matrix, write_annotation_num_matrix, write_primer_hit_evidence, write_inter_primer_dist

# Initial
from argparse import ArgumentParser
from math import floor
import sys
from Bio import SeqIO

# TODO - Go through this list and redefine/remove error messeages, and header of program. Remove the default verbose, as it is set in the argparser.
EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_INPUT_FILE_ERROR = 3
EXIT_DEPENDENCY_ERROR = 4
DEFAULT_MIN_LEN = 0
DEFAULT_VERBOSE = False
HEADER = 'FILENAME\tNUMSEQ\tTOTAL\tMIN\tAVG\tMAX'
PROGRAM_NAME = "Magphi"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="%Y-%m-%dT%H:%M:%S%z")
        logging.info('program started')
        logging.info(f"command line: {' '.join(argv)}")

        # TODO - should default logging be command-line arguements, time of run, dependency versions, Errors and progress?
        #   -log flag can be given for more info?


def main():
    start_time = time.time()

    # Retrieve the flags given by the user in the commandline
    # TODO - add in argument to not print sequences with breaks - default: Do not print
    # TODO - add in argument to not output sequences and gffs but only overview tables.
    cmd_args = get_commandline_arguments(argv[1:], PROGRAM_VERSION)


    "Orchestrate the execution of the program"
    # options = parse_args()
    init_logging(cmd_args.log)

    # Check dependencies for Magphi
    dependencies_return = check_dependencies_for_main(verbose=False) # TODO make commandline verbose controlled.
    if dependencies_return:
        logging.info("All dependencies are go!")
        logging.info("Dependency versions:")
        dependencies = ['Biopython', 'Pybedtools', 'Bedtools', 'Samtools']
        for i in range(0, len(dependencies_return)):
            logging.info(f"{dependencies[i]} v.{dependencies_return[i]}")
    else:
        warnings.warn("Some dependencies are untested versions")

    # Check the input files
    file_type, is_input_gzipped = check_inputs(cmd_args.genomes)

    # Try to construct the output folder and except if it does exist
    # TODO - make verbose controlled and log
    try:
        print("Trying to construct output folder...")
        os.mkdir(cmd_args.out_path)
        print("Succeeded!")
    except FileExistsError:
        print("Output folder exists")

    # construct a temporary folder to hold files
    # TODO - log construction of tmp folder
    tmp_folder = os.path.join(cmd_args.out_path, "Magphi_tmp_folder")
    try:
        os.mkdir(tmp_folder)
    except FileExistsError:
        raise Warning("A temporary folder already exists at the given output location. "
                      "Most likely from an incomplete analysis")

    # Check if input is GFF3 files and split genome from annotations and assign to be handed over to blast,
    # If files are not gff then assign the fastas from the input an no annotations.
    # TODO - make verbose controlled - add logging
    # TODO - should this splitting be done when each genome is being searched for primers? This will decrease the load of memory used
    if file_type == 'gff':
        print("Splitting GFF files into annotations and genomes")
        genomes, annotations = split_gff_files(cmd_args.genomes, tmp_folder, is_input_gzipped)
    else:
        genomes = cmd_args.genomes
        annotations = [None] * len(cmd_args.genomes)

    # Read in and combine primers into pairs
    primer_pairs, primer_dict = handle_primers(cmd_args.seeds)

    # Construct master dict to hold the returned information from primers
    master_primer_hits = {}
    master_annotation_hits = {}
    master_primer_evidence = {}
    primers_w_breaks = primer_pairs.copy()
    master_inter_primer_dist = {}

    # TODO - make verbose controlled - logging
    print(f'{len(genomes)} input files to be processed, starting now!')

    genomes_processed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=cmd_args.cpu) as executor:
        results = [executor.submit(screen_genome_for_primers, genomes[i], primer_pairs, cmd_args.seeds,
                                   tmp_folder, cmd_args.include_primers, file_type, annotations[i],
                                   cmd_args.out_path, cmd_args.max_primer_dist) for i, genome in enumerate(genomes)]

        for f in concurrent.futures.as_completed(results):
            genomes_processed += 1

            if genomes_processed % 25 == 0 or genomes_processed == 0:
                print(f'   File number {genomes_processed} has been processed')

            primer_hits, annots_per_interval, genome_name, primer_evidence, break_primers, inter_primer_dist = f.result()

            # Polish the genome name for the output dict:
            genome_name = genome_name.rsplit('/', 1)[-1]

            # Update the master dicts with information from current run.
            master_annotation_hits[genome_name] = annots_per_interval
            # master_annotation_hits[genome_name] = annots_per_interval
            master_primer_hits[genome_name] = primer_hits
            # master_primer_hits[genome_name] = primer_hits
            master_primer_evidence[genome_name] = primer_evidence
            master_inter_primer_dist[genome_name] = inter_primer_dist

            # Add the genome to the dict to be writen out
            master_primer_hits[genome_name]['genome'] = genome_name
            master_annotation_hits[genome_name]['genome'] = genome_name
            master_primer_evidence[genome_name]['genome'] = genome_name
            master_inter_primer_dist[genome_name]['genome'] = genome_name

            # Check if there are any primers returned as being neighbours to contig breaks.
            if len(break_primers.keys()) > 0:
                primers_w_breaks.update(break_primers)

    print('Partitioning output files into primer folders')
    # Partition output files into their primer set of origin.
    partition_outputs(primer_pairs, cmd_args.out_path)

    print('Writing output matrices')
    write_primer_hit_matrix(master_primer_hits, primer_pairs, cmd_args.out_path)
    if file_type == 'gff':
        write_annotation_num_matrix(master_annotation_hits, primers_w_breaks, cmd_args.out_path)
    write_primer_hit_evidence(master_primer_evidence, primer_pairs, cmd_args.out_path)
    write_inter_primer_dist(master_inter_primer_dist, primers_w_breaks, cmd_args.out_path)

    # Write a file specifying which primers were paired and their common name
    write_paired_primers(primer_pairs, cmd_args.out_path)

    # log quick stats for the evidence levels

    os.rmdir(tmp_folder)

    time_to_finish = time.time() - start_time
    time_to_finish = int(round(time_to_finish, 0))
    print(f"Done in: {time_to_finish} Seconds")

    # TODO - How does Log work? Do we print something to terminal?
    # TODO - Functional tests


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()

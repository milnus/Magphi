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
from math import floor, ceil
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
    from Magphi.seed_handling import handle_seeds
except ModuleNotFoundError:
    from seed_handling import handle_seeds

try:
    from Magphi.search_insertion_sites import screen_genome_for_seeds
except ModuleNotFoundError:
    from search_insertion_sites import screen_genome_for_seeds

try:
    from Magphi.wrangle_outputs import partition_outputs, write_paired_seeds
except ModuleNotFoundError:
    from wrangle_outputs import partition_outputs, write_paired_seeds

try:
    from Magphi.write_output_csv import write_seed_hit_matrix, write_annotation_num_matrix, write_seed_hit_evidence, write_inter_seed_dist
except ModuleNotFoundError:
    from write_output_csv import write_seed_hit_matrix, write_annotation_num_matrix, write_seed_hit_evidence, write_inter_seed_dist

# Initial
from argparse import ArgumentParser
from math import floor
import sys
from Bio import SeqIO

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_INPUT_FILE_ERROR = 3
EXIT_DEPENDENCY_ERROR = 4
DEFAULT_VERBOSE = False
PROGRAM_NAME = "Magphi"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def init_logging(debug_log, quiet, out_path):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        debug_log: Boolean if the logger should be a debug of info log
        out_path: Path to the output folder
    Result:
        a stream and a file logger for logging to file and stdout
    '''
    if debug_log:
        level = logging.DEBUG
    elif quiet:
        level = logging.WARNING
    else:
        level = logging.INFO

    # Construct logger logging to file
    file_logger = logging.getLogger(__name__)
    file_logger.setLevel(level)

    formatter = logging.Formatter('[%(asctime)s] %(levelname)s - %(module)s - %(message)s',
                                  datefmt="%Y-%m-%dT%H:%M:%S%z")

    file_handler = logging.FileHandler(os.path.join(out_path, 'Magphi.log'))
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    file_logger.addHandler(file_handler)

    # Log command-line argument and debug line for Magphi start
    file_logger.info(f"command line: {' '.join(argv)}")

    return file_logger


def stream_logging(file_logger):
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    file_logger.addHandler(stream_handler)

    file_logger.info('Processing started')

    return file_logger


def main():
    ''' This is the main function for running Magphi. Required arguments are genomes and seeds as a multi fasta file'''
    start_time = time.time()

    # Retrieve the flags given by the user in the commandline
    # TODO - add in argument to change tmp directory to user defined place (scratch which is faster?)
    cmd_args = get_commandline_arguments(argv[1:], PROGRAM_VERSION)

    # Try to construct the output folder and except if it does exist
    try:
        os.mkdir(cmd_args.out_path)
    except FileExistsError:
        warnings.warn("Output folder already exists")
        pass

    "Orchestrate the execution of the program"
    file_logger = init_logging(cmd_args.log, cmd_args.quiet, cmd_args.out_path)

    # Check dependencies for Magphi and logging of versions to file
    dependencies_return = check_dependencies_for_main(verbose=False)
    if dependencies_return:
        file_logger.info("Dependency versions:")
        dependencies = ['Biopython', 'Pybedtools', 'Bedtools', 'Samtools']
        for i in range(0, len(dependencies_return)):
            version = str(dependencies_return[i]).replace("\n", "")
            file_logger.info(f'{dependencies[i]} v.{version}')
    else:
        file_logger.warning("Some dependencies are untested version(s)")

    file_logger = stream_logging(file_logger)

    # Check the input files
    file_type, is_input_gzipped = check_inputs(cmd_args.genomes, file_logger)

    # construct a temporary folder to hold files
    file_logger.debug("Try to construct output folder")
    tmp_folder = os.path.join(cmd_args.out_path, "Magphi_tmp_folder")
    try:
        os.mkdir(tmp_folder)
        file_logger.debug("Output folder construction successful")
    except FileExistsError:
        file_logger.warning("A temporary folder already exists at the given output location. "
                                     "Most likely from an incomplete analysis")

    # Read in and combine seeds into pairs
    file_logger.debug("Start handling of input seed sequences")
    seed_pairs = handle_seeds(cmd_args.seeds, file_logger)

    print_seq_out = 'output' if cmd_args.no_seqs else 'None'
    if print_seq_out == 'output':
        print_seq_out = 'All' if cmd_args.print_breaks else print_seq_out

    # Construct master dict to hold the returned information from seeds
    master_seed_hits = {}
    master_annotation_hits = {}
    master_seed_evidence = {}
    seeds_w_breaks = seed_pairs.copy()
    master_inter_seed_dist = {}

    # num_genomes = len(genomes)
    num_genomes = len(cmd_args.genomes)
    # Calculate when to log progress:
    progress_num = num_genomes / 10
    progress_num = progress_num if progress_num > 1 else 1
    progress_num = ceil(progress_num)

    file_logger.info(f'{num_genomes} input files to be processed, starting now!')
    genomes_processed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=cmd_args.cpu) as executor:
        results = [executor.submit(screen_genome_for_seeds, cmd_args.genomes[i], seed_pairs, cmd_args.seeds,
                                   tmp_folder, cmd_args.include_seeds, file_type,
                                   cmd_args.out_path, cmd_args.max_seed_dist, file_logger, is_input_gzipped, print_seq_out)
                   for i, genome in enumerate(cmd_args.genomes)]

        for f in concurrent.futures.as_completed(results):
            genomes_processed += 1

            if genomes_processed % progress_num == 0 or genomes_processed == 1:
                file_logger.info(f'\tFile number {genomes_processed} has been processed')

            seed_hits, annots_per_interval, genome_name, seed_evidence, break_seeds, inter_seed_dist = f.result()

            # Polish the genome name for the output dict, and possible .fasta and .gff file extensions:
            genome_name = genome_name.rsplit('/', 1)[-1]
            genome_name = genome_name.rsplit('.gff', 1)[0]
            genome_name = genome_name.rsplit('.fasta', 1)[0]
            file_logger.debug(f'\t\tCurrently handling results from: {genome_name}')

            # Update the master dicts with information from current run.
            master_annotation_hits[genome_name] = annots_per_interval
            # master_annotation_hits[genome_name] = annots_per_interval
            master_seed_hits[genome_name] = seed_hits
            # master_seed_hits[genome_name] = seed_hits
            master_seed_evidence[genome_name] = seed_evidence
            master_inter_seed_dist[genome_name] = inter_seed_dist

            # Add the genome to the dict to be writen out
            master_seed_hits[genome_name]['genome'] = genome_name
            master_annotation_hits[genome_name]['genome'] = genome_name
            master_seed_evidence[genome_name]['genome'] = genome_name
            master_inter_seed_dist[genome_name]['genome'] = genome_name

            # Check if there are any seeds returned as being neighbours to contig breaks.
            if len(break_seeds.keys()) > 0:
                seeds_w_breaks.update(break_seeds)

    # Partition output files into their seed set of origin.
    partition_outputs(seed_pairs, cmd_args.out_path, file_logger)

    file_logger.debug('Start writing output files')
    file_logger.debug('\tStart writing hit matrix file')
    write_seed_hit_matrix(master_seed_hits, seed_pairs, cmd_args.out_path)
    if file_type == 'gff':
        file_logger.debug('\tStart writing hit extracted annotation file')
        write_annotation_num_matrix(master_annotation_hits, seeds_w_breaks, cmd_args.out_path)
    file_logger.debug('\tStart writing evidence for seed sequence pair file')
    write_seed_hit_evidence(master_seed_evidence, seed_pairs, cmd_args.out_path)
    file_logger.debug('\tStart writing distance between seed sequence pair file')
    write_inter_seed_dist(master_inter_seed_dist, seeds_w_breaks, cmd_args.out_path)
    file_logger.debug('\tStart writing seed sequence pairing file')
    write_paired_seeds(seed_pairs, cmd_args.out_path)

    # log quick stats for the evidence levels
    file_logger.debug("Remove temporary folder")
    os.rmdir(tmp_folder)

    time_to_finish = time.time() - start_time
    time_to_finish = int(round(time_to_finish, 0))
    file_logger.info(f"Magphi completed in: {time_to_finish//60} minutes {time_to_finish%60} Seconds")


if __name__ == '__main__':
    main()

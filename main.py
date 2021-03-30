import os
from sys import argv
from commandline_interface import get_commandline_arguments
from check_inputs import check_inputs
from split_gff_file import split_gff_files
from primer_handling import handle_primers
import concurrent.futures
from search_insertion_sites import screen_genome_for_primers
from write_output_csv import write_primer_hit_matrix, write_annotation_num_matrix, write_primer_hit_evidence, write_inter_primer_dist
from wrangle_outputs import partition_outputs, write_paired_primers
import time


def main():
    start_time = time.time()
    # Retrieve the flags given by the user
    cmd_args = get_commandline_arguments(argv[1:])

    # Check the types of input files
    file_type = check_inputs(cmd_args.genomes)

    # Try to construct the output folder and except if it does exist
    # TODO - make verbose controlled
    try:
        print("Trying to construct output folder...")
        os.mkdir(cmd_args.out_path)
        print("Succeeded!")
    except FileExistsError:
        print("Output folder exists")

    # construct a temporary folder to hold files
    tmp_folder = os.path.join(cmd_args.out_path, "phupa_tmp_folder")
    try:
        os.mkdir(tmp_folder)
    except FileExistsError:
        raise Warning("A temporary folder already exists at the given output location. "
                      "Most likely from an incomplete analysis")

    # Check if input is GFF3 files and split genome from annotations and assign to be handed over to blast,
    # If files are not gff then assign the fastas from the input an no annotations.
    if file_type == 'gff':
        # TODO - Possibly parallellise this process
        genomes, annotations = split_gff_files(cmd_args.genomes, tmp_folder)
    else:
        genomes = cmd_args.genomes
        annotations = [None]*len(cmd_args.genomes)

    # Read in and combine primers into pairs
    primer_pairs, primer_dict = handle_primers(cmd_args.primers)

    """
    OLD SINGLE THREAD CODE THAT HAS BEEN SUBSTITUTED WITH PARALLEL CODE BELOW
    # for i, genome in enumerate(genomes):
    #     primer_hits, annots_per_interval = screen_genome_for_primers(genomes[i], primer_pairs, cmd_args.primers, tmp_folder, cmd_args.include_primers,
    #                               file_type, annotations[i], cmd_args.out_path)
    #
    #     print(primer_hits)
    #     print(annots_per_interval)
    """

    # Construct master dict to hold the returned information from primers
    master_primer_hits = {}
    master_annotation_hits = {}
    master_primer_evidence = {}
    primers_w_breaks = primer_pairs.copy()
    master_inter_primer_dist = {}

    # TODO - make verbose controlled
    print(f'{len(genomes)} input files has to be processed, starting now!')

    # TODO - add a command line argument that allows the user to set the number of workers.
    with concurrent.futures.ProcessPoolExecutor(max_workers=cmd_args.cpu) as executor:
        results = [executor.submit(screen_genome_for_primers, genomes[i], primer_pairs, cmd_args.primers,
                                   tmp_folder, cmd_args.include_primers, file_type, annotations[i],
                                   cmd_args.out_path, cmd_args.max_primer_dist, i) for i, genome in enumerate(genomes)]

        for f in concurrent.futures.as_completed(results):
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

    # Partition output files into their primer set of origin.
    partition_outputs(primer_pairs, cmd_args.out_path)

    # TODO - Output the Master dicts as matrices.
    #   * Input the primer names to be used with replace to remove them from the genome name and give columns the primers as name.
    #   * Possibly output the length of sequences found between primers along with the other information.
    # write_output_matrixes(master_primer_hits, master_annotation_hits, primer_pairs, cmd_args.out_path)
    write_primer_hit_matrix(master_primer_hits, primer_pairs, cmd_args.out_path)
    if file_type == 'gff':
        write_annotation_num_matrix(master_annotation_hits, primers_w_breaks, cmd_args.out_path)
    write_primer_hit_evidence(master_primer_evidence, primer_pairs, cmd_args.out_path)
    write_inter_primer_dist(master_inter_primer_dist, primers_w_breaks, cmd_args.out_path)

    # Write a file specifying which primers were paired and their common name
    write_paired_primers(primer_pairs, cmd_args.out_path)


    # TODO - Format output in ways that are usefull. Prokka style gff3 format with annotations and sequnce, keep them seperate? genbank? maybe only print the annotation output if there is a gene in the region

    os.rmdir(tmp_folder)

    time_to_finish = time.time() - start_time
    time_to_finish = int(round(time_to_finish, 0))
    print(f"Done in: {time_to_finish} Seconds")

if __name__ == '__main__':
    main()

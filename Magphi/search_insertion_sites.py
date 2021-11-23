import os
import glob
import csv
import warnings
import sys
import gzip
from itertools import combinations
from shutil import copyfile, copyfileobj
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Sequencing.Applications import SamtoolsFaidxCommandline
from Bio import SearchIO
import pybedtools as bedtools

try:
    from Magphi.split_gff_file import split_single_gff
except ModuleNotFoundError:
    from split_gff_file import split_single_gff
# pylint: disable=E1123


def blast_insertion_site(seeds, genome_file, tmp_name):
    # Construct genome_db path and name
    genome_db = f'{genome_file}_tmp_db'
    # Make blast database command for given genome
    c_line_makedb = NcbimakeblastdbCommandline(dbtype='nucl', input_file=genome_file, out=genome_db)

    # Run makeblastdb in command line
    c_line_makedb()

    blast_out_xml = f'{tmp_name}.xml'
    c_line = NcbiblastnCommandline(query=seeds,
                                   db=genome_db,
                                   evalue=0.001,
                                   outfmt=5,
                                   qcov_hsp_perc=50,
                                   out=blast_out_xml)

    # Run the blast command in commandline blast
    c_line()

    # Delete blast database
    file_list = glob.glob(f'{genome_file}_tmp_db.*')
    for file in file_list:
        os.remove(file)

    return blast_out_xml


def blast_out_to_sorted_bed(blast_xml_output, include_seeds, genome_name, seed_pairs):
    # Read in the blast output into a generator object
    blast_output = SearchIO.parse(blast_xml_output, format='blast-xml')

    # Construct list to hold lines for seed bed files
    bed_list = []
    for qresult in blast_output:
        bed_line = [None]*4
        # Insert the name of the seed
        bed_line[3] = qresult.id
        for hit in qresult:
            bed_line_hit = bed_line.copy()
            # Get the contig that has been hit by the seed
            bed_line_hit[0] = hit.id
            for hsp in hit:
                bed_line_hsp = bed_line_hit.copy()
                # Extract the beginning of the seed (extract on to make 0-based for BED file)
                bed_line_hsp[1] = hsp.hit_start
                # Extract the end of the seed (extract on to make 0-based for BED file)
                bed_line_hsp[2] = hsp.hit_end

                # Add the bed line for the hit to the bed file list
                bed_list.append(bed_line_hsp)

    # Construct lists to hold bed file names, and a dict to hold the number of hits for each seed pair
    blast_hit_beds = []
    exclusion_list = []
    seed_hits = {}

    # Split up the bed lines into each seed and write a bed file for each seed
    for seed_pair in seed_pairs.keys():
        # Find bed lines for seed pair
        seed_bed_lines = [line for line in bed_list
                            if seed_pairs[seed_pair][0] == line[3] or seed_pairs[seed_pair][1] == line[3]]

        # Record the number of hits in a genome by a set of seeds
        seed_hits[seed_pair] = len(seed_bed_lines)

        # Write the bed files for each seed pair
        bed_name = f'{genome_name}~~{seed_pair}.bed'
        with open(bed_name, 'w') as bed_file:
            writer = csv.writer(bed_file, delimiter='\t')
            for line in seed_bed_lines:
                writer.writerow(line)
        bed_file.close()
        blast_hit_beds.append(bed_name)

        # Sort the bed file
        sorted_bed = bedtools.BedTool(bed_name).sort()
        sorted_bed.saveas(bed_name)

        # Check if seeds are to be excluded, if then copy seed bed to be used later in bedtools intersect
        if not include_seeds:
            copyfile(f'{genome_name}~~{seed_pair}.bed', f'{genome_name}~~{seed_pair}_seeds.bed')
            exclusion_list.append(f'{genome_name}~~{seed_pair}_seeds.bed')

    return blast_hit_beds, exclusion_list, seed_hits


def write_bed_from_list_of_seeds(list_of_seeds, bed_file_name):
    # Write out new .bed file with only the extended seeds.
    # Join the items in the lists (seeds) into a sting
    list_of_seeds = [' '.join(seed) for seed in list_of_seeds]

    # Join the seeds into one large sting
    list_of_seeds = '\n'.join(list_of_seeds)

    # Construct BedTools object
    bed_of_seeds = bedtools.BedTool(list_of_seeds, from_string=True)

    # Sort the BedTools object
    bed_of_seeds = bed_of_seeds.sort()

    # Save the seeds into a .bed file
    bed_of_seeds.saveas(bed_file_name)


def seed_reach_contig_end_calc(genome_file, max_seed_dist, seed_contig_hits):
    # Construct dict to hold the indexing of the fasta file
    fai_dict = {}
    # Read and search the .fai file to get the edge of contigs
    with open(genome_file, 'r') as fai_file:
        for line in fai_file.readlines():
            line = line.split('\t')
            fai_dict[line[0]] = int(line[1])
    fai_file.close()

    # Extract all lines from dict holding Bed intervals for seed hits
    intervals = [interval for key in seed_contig_hits.keys() for interval in seed_contig_hits[key]]

    # Insert line identifier for interaction matrix
    [line.append(str(i)) for i, line in enumerate(intervals)]

    # Construct matrix to hold info on contig ends reached
    end_reached_matrix = [[0] * len(intervals) for _ in range(len(intervals))]

    # Go through each seed hit to see if it reached the end of the contig under the max distance
    for i, interval in enumerate(intervals):
        # Check if the 5' and 3' can be reached using the max distance
        if int(interval[1]) - max_seed_dist < 0:
            end_reached_matrix[i][0] = 1

        if int(interval[2]) + max_seed_dist > fai_dict[interval[0]]:
            end_reached_matrix[i][1] = 1

    # Find how many ends each seed can reach
    end_sums = [sum(ends) for ends in end_reached_matrix]

    # Count the number of seeds that reach one or both ends
    end_reaches = [end_sums.count(2), end_sums.count(1)]

    return end_reaches, end_sums, end_reached_matrix, intervals


def examine_flanking_regions(seed_contig_hits, max_seed_dist, genome_file, bed_file_name=None):
    '''
    Function that takes a set of seeds on contigs, and search for how they can be connected using a genome and length of contigs.
    :param seed_contig_hits: dict with each key being a contig and values being a list of seed sequence hit on contigs as lists
    :param max_seed_dist: int for the maximum distance allowed between seed sequenes
    :param genome_file: .fai (fasta-index) of a genome used as info on contig sizes
    :param bed_file_name: Used for saving updated BED file with coordinates
    :return: Evidence level for a set of seed seuqences that have hit the genome.
    '''
    # If no limits are given and only one contig is hit then no pairs can be predicted
    if max_seed_dist == 0:
        return 2

    # If limits are given then check overlap of seeds, if they are only found on one contig
    elif len(seed_contig_hits) == 1 and max_seed_dist > 0:
        # Extract the intervals found for seeds
        intervals = seed_contig_hits[list(seed_contig_hits.keys())[0]]

        # Construct an interaction matrix for the seeds
        interaction_matrix = [[0]*len(intervals) for _ in range(len(intervals))]

        # Insert line identifier for interaction matrix
        [line.append(str(i)) for i, line in enumerate(intervals)]

        for i, interval in enumerate(intervals):
            # temp intervals to be altered
            intervals_tmp = intervals.copy()

            # Convert interval from list to BedTools object
            bed_interval = ' '.join(interval)
            bed_interval = bedtools.BedTool(bed_interval, from_string=True)

            # slop on the object with the max distance on both sides
            bed_interval = bed_interval.slop(b=max_seed_dist, g=genome_file)

            # Reinsert the newly slopped interval
            intervals_tmp[i] = list(bed_interval[0])

            # Construct BedTool object from all seed hits in intervals
            # Join the intervals into strings
            intervals_tmp = [' '.join(hit) for hit in intervals_tmp]
            # Join the string intervals into one big string
            intervals_tmp = '\n'.join(intervals_tmp)

            # Create a BedTool object from the large string
            intervals_bed = bedtools.BedTool(intervals_tmp, from_string=True)

            # Resort if needed
            intervals_bed = intervals_bed.sort()

            # merge
            merged_intervals = intervals_bed.merge(d=0, c=[4, 4, 5], o='count,collapse,collapse')

            # Check to see if some lines have been merged, if then record this in the interaction matrix
            interacting_lines = [list(line) for line in merged_intervals if int(line[3]) >= 2]

            # Filter out overlapping seeds in which the seeds are the same
            interacting_ids = [line[5] for line in interacting_lines if len(set(line[4].split(','))) >= 2]
            interacting_seed_names = [line[4] for line in interacting_lines if len(set(line[4].split(','))) >= 2]

            # Insert interacting seeds in matrix (only in upper part of matrix):
            if len(interacting_ids) > 0:
                for cluster in zip(interacting_ids, interacting_seed_names):
                    # split the pairs in the zipped tuple and find the combinations of names and ids
                    id_combinations = combinations(cluster[0].split(','), 2)
                    name_combinations = combinations(cluster[1].split(','), 2)

                    # Go through each pair of id and name combination to record valid hits
                    for pair in zip(id_combinations, name_combinations):
                        # Check that the combination of seeds are not the same seed and that seed being examined is involved
                        if pair[1][0] != pair[1][1] and str(i) in pair[0]:
                            interaction_matrix[int(pair[0][0])][int(pair[0][1])] = 1


        # Summarise all rows of the interaction matrix
        num_interactions = [sum(row) for row in interaction_matrix]
        # Summarise col values to get a number of total interactions
        num_interactions = sum(num_interactions)

        if num_interactions == 0:
            return 1

        elif num_interactions == 1:
            # Check if seeds reach end of contigs
            end_reaches, *_ = seed_reach_contig_end_calc(genome_file, max_seed_dist, seed_contig_hits)
            if sum(end_reaches) > 0:
                return 2

            interaction_indexes = []
            # Get index of the interaction between seeds
            for i, interacts in enumerate(interaction_matrix):
                if 1 in interacts:
                    for j, interaction in enumerate(interaction_matrix[i]):
                        if interaction == 1:
                            interaction_indexes = [i, j]

            # Extract the seeds that correspond to the identified indexes
            interaction_seeds = [intervals[index][:4] for index in interaction_indexes]
            # Write the bed file from the list of seeds
            write_bed_from_list_of_seeds(interaction_seeds, bed_file_name)
            return '5B'

        elif num_interactions > 1:
            return 2

    else:
        end_reaches, end_sums, end_reached_matrix, intervals = seed_reach_contig_end_calc(genome_file, max_seed_dist, seed_contig_hits)

        # See if no ends has been reached
        if end_reaches[0] == 0 and end_reaches[1] == 0:
            return 1

        # Check if one seed has reached one or both ends, while the other has reached no ends.
        elif (end_reaches[0] == 1 and end_reaches[1] == 0) or (end_reaches[0] == 0 and end_reaches[1] == 1):
            return 3

        # Check if
        # two or more seeds hit both ends of their contig (end_reaches[0] > 1) or
        # three or more seeds all hit one end of their contig (end_reaches[1] > 2) or
        # one seed reaches two ends and two seeds reaches one end each (end_reaches[0] == 1 and end_reaches[1] == 2)
        elif end_reaches[0] > 1 or end_reaches[1] > 2 or (end_reaches[0] == 1 and end_reaches[1] == 2):
            return 2

        # Check if one seed has reached both ends and one seed has reached one end
        elif end_reaches[1] == 1 and end_reaches[0] == 1:
            return 2

        # Check if two seeds have reached one end.
        elif end_reaches[1] == 2:
            # Organise the bed and reprint with the usefull seeds and adjusted coordinates to end
            seed_indexes = [index for index, reached in enumerate(end_sums) if reached == 1]

            # Initialise list to hold extended seeds
            extended_seeds = []

            for index in seed_indexes:
                # Get the seed extension hits
                seed_extension_hits = end_reached_matrix[index]

                # Get the seed in question, and remove the identifier at the end.
                extend_seed = intervals[index][:4]

                # Convert seed to a BedTools object
                extend_seed = bedtools.BedTool(' '.join(extend_seed), from_string=True)

                # Examine if seed should be extended in 5' end
                if seed_extension_hits[0] == 1:
                    extend_seed = extend_seed.slop(g=genome_file, l=max_seed_dist, r=0)

                # Examine if the seed should be extended in the 3' end
                elif seed_extension_hits[1] == 1:
                    extend_seed = extend_seed.slop(g=genome_file, r=max_seed_dist, l=0)

                else:
                    warnings.warn("A seed found to reach the end of a contig's "
                                  "end was not extended in any direction!\n"
                                  f"Please report this along with the genome: {genome_file},"
                                  f" the seed in question: {extend_seed}, "
                                  f"and the max distance set: {max_seed_dist}")
                    return "Something unexpected was found, please check warning log!"

                extended_seeds.append(extend_seed[0])

            # Write out new .bed file with only the extended seeds.
            write_bed_from_list_of_seeds(extended_seeds, bed_file_name)
            return '4B'

        else:
            warnings.warn(f'Some unaccounted for constalation of seeds hits '
                          f'was found to reach the end of contigs.\n'
                          f'Genome .fai file examined: {genome_file}.\n'
                          f'Please attach the .fai file (found in the tmp folder)'
                          f' when reporting this, along with the following: '
                          f'{intervals} and the max distance used: {max_seed_dist}')
            return "Something unexpected was found, please check warning log!"


def check_seeds_placement(bed_files, seed_pairs, seed_hits, max_seed_dist, genome_file, file_type, tmp_folder):
    ''' Function to determine the placement of seeds and how to search for connections between them
    Returns a genome file in the temporary directory and an evidence score for each bed file and its seed connections'''
    # Produce genome index file using samtools
    samtools_faidx_cmd = SamtoolsFaidxCommandline(ref=genome_file)
    samtools_faidx_cmd()

    # Initialise the dict to hold the support for a seed hit
    seed_hit_support_dict = dict.fromkeys(seed_pairs)

    # Check each BED file
    for file in bed_files.copy():
        # Strip file name to get name of genome
        seed_name = file.rsplit('.', 1)[0]
        seed_name = seed_name.rsplit('~~', 1)[-1]

        # Initialise dict to hold information on seed distribution over contigs
        seed_to_contig = {}

        # Check if more than one seeds has hit the genome,
        # if, then check if they are same contig
        if seed_hits[seed_name] > 1:
            # Read Bed file
            bed_file = bedtools.BedTool(file)

            # Index lines in BED file and which contig seed sequences hit
            for line in bed_file:
                contig = line[0]
                try:
                    # Add the line to the contig
                    seed_to_contig[contig] += [list(line)]
                except KeyError:
                    # If not found add the contig and the line
                    seed_to_contig[contig] = [list(line)]

            # Find number of different seeds that hit the contigs
            # Get all hits on contigs an unpack them into a single list of lists, then count the number of unique seed sequence hits
            seed_bed_lines = [seed_to_contig[key] for key in seed_to_contig.keys()]
            seed_bed_lines = [hit for contig in seed_bed_lines for hit in contig]
            uniq_seeds = len(set([hit[3] for hit in seed_bed_lines]))

            # If more than one contig is hit and more than one seed sequence has hit then examine further,
            # else only one contig is hit
            if len(seed_to_contig.keys()) > 1 and uniq_seeds > 1:
                # Examine if some of the seeds are on the same contig:
                multi_hit_contigs = [contig for contig in seed_to_contig.keys() if len(seed_to_contig[contig]) >= 2]
                # Check if there are any contigs hit with multiple seeds,
                # if, then examine them for overlaps between seeds,
                # if not, examine across contigs
                if len(multi_hit_contigs) > 0:
                    # Examine all contigs that contain multiple seeds for an overlap
                    for contig in multi_hit_contigs:
                        hit_contig = {contig: seed_to_contig[contig]}

                        return_value = examine_flanking_regions(hit_contig, max_seed_dist, f'{genome_file}.fai', file+'_5')

                        # Record the return value as evidence
                        # 1 or 2 means more examination is required
                        if return_value == 2:
                            seed_hit_support_dict[seed_name] = return_value
                            bed_files.remove(file)
                        elif return_value == 1:
                            # Examine the remaining junctions to see if two seeds can be found to connect across contigs
                            seed_hit_support_dict[seed_name] = examine_flanking_regions(seed_to_contig,
                                                                                            max_seed_dist,
                                                                                            f'{genome_file}.fai',
                                                                                            file)
                            # Check if no connection could be made on or across contigs, if then delete bed from further processing
                            if seed_hit_support_dict[seed_name] == 1 or seed_hit_support_dict[seed_name] == 2:
                                bed_files.remove(file)

                        elif return_value == '5B':
                            seed_hit_support_dict[seed_name] = return_value
                            alternative_return_value = examine_flanking_regions(seed_to_contig,
                                                                                max_seed_dist,
                                                                                f'{genome_file}.fai',
                                                                                file+'_6')

                            if alternative_return_value == '4B': # ADD alternative value 3?
                                seed_hit_support_dict[seed_name] = 2
                                os.remove(file+'_5')
                                try:
                                    os.remove(file + '_6')
                                except FileNotFoundError:
                                    pass

                            else:
                                seed_hit_support_dict[seed_name] = return_value
                                os.rename(file+'_5', file)

                        else:
                            sys.exit(100)
                else:
                    # Examine the seed hits across different contigs
                    return_value = examine_flanking_regions(seed_to_contig,
                                                            max_seed_dist,
                                                            f'{genome_file}.fai',
                                                            file)
                    if return_value == 1 and seed_hits[seed_name] == 2 and uniq_seeds == 2:
                        seed_hit_support_dict[seed_name] = '4A'
                    elif return_value == 2:
                        seed_hit_support_dict[seed_name] = return_value
                        bed_files.remove(file)
                    else:
                        seed_hit_support_dict[seed_name] = return_value

            else:
                # Check if precisely two seeds hit the contig and that the mates from the pair hit one time each
                if seed_hits[seed_name] == 2 and uniq_seeds == 2:
                    # Test if seeds can reach end of contig
                    return_value = examine_flanking_regions(seed_to_contig, max_seed_dist, f'{genome_file}.fai')
                    # Score seed hit
                    if return_value == 1:
                        seed_hit_support_dict[seed_name] = '5A'
                    else:
                        seed_hit_support_dict[seed_name] = return_value

                # Check if two unique seeds hit, but one/both may have hit multiple times
                elif seed_hits[seed_name] > 2 and uniq_seeds == 2:
                    # Examine seed hits and record returned score
                    seed_hit_support_dict[seed_name] = examine_flanking_regions(seed_to_contig,
                                                                                    max_seed_dist,
                                                                                    f'{genome_file}.fai')
                    if seed_hit_support_dict[seed_name] == 1 or seed_hit_support_dict[seed_name] == 2:
                        bed_files.remove(file)

                # Check that only one seed has hit, but it has hit multiple times
                elif uniq_seeds == 1:
                    # Score seed hit
                    seed_hit_support_dict[seed_name] = 0
                else:
                    raise NotImplementedError(f'Some unaccounted for constalation of seeds hits '
                                              f'was found to hit a single contig.\n'
                                              f'The seed pair in question is {seed_name} in genome {genome_file}.\n'
                                              f'Please report this along with the following: {seed_bed_lines}')
        else:
            # set seed evidence and remove bed-file from further processing
            seed_hit_support_dict[seed_name] = 0
            bed_files.remove(file)

    return genome_file, seed_hit_support_dict


def bed_merge_handling(blast_hit_beds, include_seeds, exclude_seed_list, max_seed_dist, seed_evidence):
    # initialise list to hold merged bed file names:
    merged_bed_files = []

    # Process the bed file for each seed pair
    for i, bed_file in enumerate(blast_hit_beds):
        seed_hits = bedtools.BedTool(bed_file)

        # Get name of seed
        split_bed_name = bed_file.rsplit('/', 1)[-1]
        split_bed_name = split_bed_name.replace('.bed', '')
        seed_name = split_bed_name.rsplit('~~', 1)[-1]

        # Check if max distance between seeds is set to unlimited (0)
        if max_seed_dist == 0:
            max_seed_dist = 9999999999

        # Merge the bed file and collapse the column containing the names and count the number of intervals collapsed
        if len(seed_hits) > 0:
            seed_hits = seed_hits.merge(c=[4, 4], o='collapse,count', d=max_seed_dist)
        else:
            seed_evidence[seed_name] = 0
            continue

        # Check if the evidence level for the seeds should increase due to seeds being merged.
        # Check that there is one line and that two seeds have been merged and not just one seed hit
        if len(seed_hits) == 1 and int(seed_hits[0][4]) == 2:
            seed_evidence[seed_name] = '5B'

        # Evaluate if seeds are to be excluded from the intervals
        if not include_seeds:
            # Load bed containing seeds to be excluded
            exclusion_bed = bedtools.BedTool(exclude_seed_list[i])

            pre_deletion_intervals = len(seed_hits)

            # Remove the seed intervals
            seed_hits = seed_hits.subtract(exclusion_bed)

            if pre_deletion_intervals > len(seed_hits) and (str(seed_evidence[seed_name]) < '3' or seed_evidence[seed_name] == '5B' or seed_evidence[seed_name] == '4B'):
                seed_evidence[seed_name] = 3

        # Save the merged intervals
        # Construct name
        merged_bed_file = bed_file.rsplit('.', 1)[0]
        merged_bed_file = f'{merged_bed_file}_merged.bed'
        # Save
        # Check if there is any interval to report,
        # else notify by setting the return evidence level
        if len(seed_hits):
            seed_hits.saveas(merged_bed_file)
            # Add to list
            merged_bed_files.append(merged_bed_file)
        else:
            if str(seed_evidence[seed_name]) < '3':
                seed_evidence[seed_name] = 3

    return merged_bed_files, seed_evidence


def extract_seqs_n_annots(merged_bed_files, file_type, genome_file, annotation_file, tmp_folder, out_path, seed_pairs,
                          seed_evidence, print_seq_out):
    """"""
    # Initiate dict to hold number of annotations per interval
    annots_per_interval = dict.fromkeys(seed_pairs)
    # Initiate dict to hold seeds that may neighbour a sequence break
    break_seeds = {}
    # Initiate dict tot hold distances between seeds
    inter_seed_dist = {}

    # Extract sequence from Fasta from all merged bed files
    for merged_bed in merged_bed_files:
        # Read Bed file
        merged_intervals = bedtools.BedTool(merged_bed)

        # Check if the bed file is empty indicate that no intervals are present by giving nan in annotations
        # There are no intervals when excluding seeds and no seeds are merged for specific seed set.
        if file_type == 'gff':
            seed_pair = merged_bed.rsplit('~~', 1)[-1]
            seed_pair = seed_pair.rsplit('_merged', 1)[0]
            if len(merged_intervals) == 0:
                annots_per_interval[seed_pair] = 'nan'
                continue

        # Go though each line of bed, really only necessary if there are breaks in genome.
        for interval in merged_intervals:
            # Identify seed pair for interval:
            seeds = interval[3].split(',')
            # Check if seeds are merged, or if single seed interval
            if len(seeds) == 2:
                seed_pair_name = [seed_name for seed_name in seed_pairs if
                                    seeds[0] in seed_pairs[seed_name] and seeds[1] in seed_pairs[seed_name]]

                seed_pair_name = seed_pair_name[0]
            else:
                seed_pair_name = f'{seeds[0]}_break'
                break_seeds[seed_pair_name] = [seeds[0], 'break']

            # Construct output name for extracted sequence
            genome_name = genome_file.rsplit('/', 1)[-1]
            genome_name = genome_name.rsplit('.', 1)[0]
            genome_name = genome_name.rsplit('_tmp', 1)[0]
            out_file_name = f'{genome_name}--{seed_pair_name}'
            output_genome = os.path.join(out_path, f'{out_file_name}.fasta')

            # Construct good fasta header
            fasta_header = f'{genome_name}_{seed_pair_name}'

            # Convert interval to string to be made into a BedTool object
            interval = f'{interval[0]} {int(interval[1])} {int(interval[2])} {fasta_header}'

            # Create BedTool object from line
            interval = bedtools.BedTool(interval, from_string=True)

            # Record the length of the segments between seeds
            inter_seed_dist[seed_pair_name] = len(interval[0])

            # Extract the fasta sequence and save it in output folder
            if print_seq_out == 'All' or '_break' not in seed_pair_name and print_seq_out == 'output':
                try:
                    interval.sequence(fi=genome_file, fo=output_genome, nameOnly=True)
                except bedtools.helpers.BEDToolsError:
                    interval.sequence(fi=genome_file, fo=output_genome, name=True)

            # extract annotations if gff is provided as input
            if file_type == 'gff':
                # Load in the gff as a BedTool object
                gff_file = bedtools.BedTool(annotation_file)

                # Save the merged BED interval currently being handled
                interval_save = os.path.join(tmp_folder, f'interval_{genome_name}')
                interval.saveas(interval_save)

                # extract the annotations within the current merged bed interval
                annot = gff_file.intersect(b=interval, f=0.95)

                # Record the number of annotations within each interval
                annots_per_interval[seed_pair_name] = len(annot)

                # Start writing the output file, if an annotation if found
                if len(annot) > 0:
                    # Check if seed sequence reach a contig break
                    # If then get the seed pair to adjust the evidence level
                    if '_break' in seed_pair_name:
                        seed_pair_key = seed_pair_name
                        seed_pair_key = seed_pair_key.rsplit('_', 1)[0]
                        seed_pair_key = [key for key in seed_pairs.keys() if seed_pair_key in seed_pairs[key]][0]
                    else:
                        seed_pair_key = seed_pair_name

                    # Increase the evidence level for the seed sequence pair.
                    # Only if the evidence level is 4B or 5B and seed pair is valid (no breaks)
                    if seed_pair_key in seed_evidence.keys():
                        if seed_evidence[seed_pair_key] == '4B':
                            seed_evidence[seed_pair_key] = '4C'
                        elif seed_evidence[seed_pair_key] == '5B':
                            seed_evidence[seed_pair_key] = '5C'

                    # Construct name for output file
                    output_gff = os.path.join(out_path, f'{out_file_name}.gff')

                    # Find start coordinate for current merged interval
                    start_coordinate = int(interval[0][1])+1

                    # Write GFF output, if allowed or break
                    if print_seq_out == 'All' or '_break' not in seed_pair_name and print_seq_out == 'output':
                        # Open output file
                        gff_output_file = open(output_gff, 'w')
                        gff_writer = csv.writer(gff_output_file, delimiter='\t')

                        # Insert the required '##gff-version 3' line
                        gff_output_file.write('##gff-version 3\n')
                        # Write line recording the length of the extracted region
                        gff_output_file.write(f'##sequence-region {fasta_header} 1 {len(interval[0])}\n')

                        # Go through each line adjust the coordinates, the contig name and write the gff
                        for i, line in enumerate(annot):
                            line[0] = fasta_header
                            line[3] = int(line[3]) - start_coordinate
                            line[4] = int(line[4]) - start_coordinate + 1
                            gff_writer.writerow(line)

                        # Add in the tag for the Fasta
                        gff_output_file.write('##FASTA\n')
                        # Open the fasta sequence for copying
                        fasta_file = open(output_genome, 'r')
                        # Append the fasta sequence for the extracted interval
                        copyfileobj(fasta_file, gff_output_file)

                        # Close output files
                        gff_output_file.close()
                        fasta_file.close()

                # Delete the file containing the current merged bed interval
                os.remove(interval_save)

    # Remove the temporary genome and the index file (.fai), or the temporary fasta and annotation file for gff input
    if file_type == 'fasta':
        os.remove(genome_file)
        try:
            os.remove(f'{genome_file}.fai')
        except FileNotFoundError:
            pass
    else:
        try:
            os.remove(f'{genome_file}.fai')
        except FileNotFoundError:
            pass
        os.remove(genome_file)
        os.remove(annotation_file)

    return annots_per_interval, break_seeds, seed_evidence, inter_seed_dist


def screen_genome_for_seeds(genome_file, seed_pairs, seed_path, tmp_folder,
                              include_seeds, file_type, out_path, max_seed_dist, file_logger,
                              is_input_gzipped, print_seq_out):
    """
    Function that summarise the search of seeds sequences, determination of position and extraction.
    :param genome_file: Path to genome to be searched for seed sequences
    :param seed_pairs: Dict of seed sequences that are matched into pairs
    :param seed_path: Path to the file containing the seed sequence fasta
    :param tmp_folder: The working temporary folder Magphi has created
    :param include_seeds: Boolean to indicate if the coordinates of the seed sequences should be included in results
    :param file_type: Variable indicating if input genomes are gff or fasta files
    :param annotation_file: Path to the annotations part of an input gff file, previously split from the genome part
    :param out_path: Path to the output folder
    :param max_seed_dist: Integer indicating the maximum distance allowed between two seed sequences
    :param file_logger: File to with debug log statements should be written.
    :param is_input_gzipped: Bool to tell if the input file is gzipped
    :param print_seq_out: Bool to indicate if outputs related to fasta and gff sequences should be given
    :return: Number of times a seed sequence pairs hit the genome,
    number of annotations in the interval found between seed sequences,
    name of the genome extracted from,
    evidence for how seed sequences hit and connected in the genome,
    Any seed sequences that were next to a sequences break
    distance between seed sequences that could be connected.
    """
    file_logger.debug(f"Start search of sequences is {genome_file}")

    genome_name = os.path.basename(genome_file)
    genome_name = genome_name.rsplit('.', 1)[0]

    tmp_genome_folder = os.path.join(tmp_folder, genome_name)
    os.mkdir(tmp_genome_folder)

    if file_type == 'fasta':
        tmp_genome = os.path.join(tmp_genome_folder, genome_file.rsplit('/')[-1])

        if is_input_gzipped:
            file_logger.debug("\tCopying gzipped fasta file into tmp dir")
            tmp_genome = tmp_genome.split('.gz')[0]
            with gzip.open(genome_file, 'rt') as in_file:
                with open(tmp_genome, 'w') as out_file:
                    for line in in_file:
                        out_file.write(line)
        else:
            file_logger.debug("\tCopying ungzipped fasta file into tmp dir")
            copyfile(genome_file, tmp_genome)

        genome_file = tmp_genome
        annotation_file = None

    if file_type == 'gff':
        genome_file, annotation_file = split_single_gff(genome_file, tmp_genome_folder, is_input_gzipped)

    # Concatenate the genome name and the path to the temporary folder
    # genome_name = os.path.join(tmp_folder, genome_name)
    genome_name = os.path.join(tmp_genome_folder, genome_name)

    # Run blast with genome and insertion site sequences
    file_logger.debug(f"\tBLASTing: {genome_file}")
    blast_xml_output = blast_insertion_site(seed_path, genome_file, f'{genome_name}_blast')

    # Construct a bedfile from blast output
    file_logger.debug(f"\tConstructing bed file: {genome_file}")
    blast_hit_beds, exclude_seed_list, seed_hits = blast_out_to_sorted_bed(blast_xml_output,
                                                                               include_seeds,
                                                                               genome_name,
                                                                               seed_pairs)

    # Examine the seed hits and try to find solutions when multiple seeds hit at once
    file_logger.debug(f"\tChecking seed sequence placement: {genome_file}")
    genome_file, seed_evidence = check_seeds_placement(blast_hit_beds, seed_pairs, seed_hits, max_seed_dist,
                                                           genome_file, file_type, tmp_folder)

    # sort and merge the bed files
    file_logger.debug(f"\tMerging seed sequences: {genome_file}")
    merged_bed_files, seed_evidence = bed_merge_handling(blast_hit_beds, include_seeds, exclude_seed_list,
                                                           max_seed_dist, seed_evidence)

    # Extract sequences and annotations using merged intervals.
    file_logger.debug(f"\tExtracting sequences from intervals: {genome_file}")
    annots_per_interval, break_seeds, seed_evidence, inter_seed_dist = extract_seqs_n_annots(merged_bed_files,
                                                                                                   file_type,
                                                                                                   genome_file,
                                                                                                   annotation_file,
                                                                                                   tmp_folder,
                                                                                                   out_path,
                                                                                                   seed_pairs,
                                                                                                   seed_evidence,
                                                                                                   print_seq_out)

    # clean up by removing blast xml output, and merged and seed hit beds from tmp folder:
    os.remove(blast_xml_output)
    tmp_folder_files = os.listdir(tmp_genome_folder)
    for file in tmp_folder_files:
        if str(genome_name).rsplit('/', 1)[1] in file:
            try:
                os.remove(os.path.join(tmp_genome_folder, file))
            except FileNotFoundError:
                pass
    os.rmdir(tmp_genome_folder)
    # tmp_folder_files = os.listdir(tmp_folder)
    # for file in tmp_folder_files:
    #     if '.bed' in file and str(genome_name).rsplit('/')[1] in file:
    #         os.remove(os.path.join(tmp_folder, file))

    return seed_hits, annots_per_interval, genome_name, seed_evidence, break_seeds, inter_seed_dist

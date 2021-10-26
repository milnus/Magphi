import os
import glob
import csv
import warnings
import sys
from itertools import combinations
from shutil import copyfile, copyfileobj
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Sequencing.Applications import SamtoolsFaidxCommandline
from Bio import SearchIO
import pybedtools as bedtools
# pylint: disable=E1123


def blast_insertion_site(primers, genome_file, tmp_name):
    # Construct genome_db path and name
    genome_db = f'{genome_file}_tmp_db'
    # Make blast database command for given genome
    c_line_makedb = NcbimakeblastdbCommandline(dbtype='nucl', input_file=genome_file, out=genome_db)

    # Run makeblastdb in command line
    c_line_makedb()

    blast_out_xml = f'{tmp_name}.xml'
    c_line = NcbiblastnCommandline(query=primers,
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


def blast_out_to_sorted_bed(blast_xml_output, include_primers, genome_name, primer_pairs):
    # Read in the blast output into a generator object
    blast_output = SearchIO.parse(blast_xml_output, format='blast-xml')

    # Construct list to hold lines for primer bed files
    bed_list = []
    for qresult in blast_output:
        bed_line = [None]*4
        # Insert the name of the primer
        bed_line[3] = qresult.id
        for hit in qresult:
            bed_line_hit = bed_line.copy()
            # Get the contig that has been hit by the primer
            bed_line_hit[0] = hit.id
            for hsp in hit:
                bed_line_hsp = bed_line_hit.copy()
                # Extract the beginning of the primer (extract on to make 0-based for BED file)
                bed_line_hsp[1] = hsp.hit_start
                # Extract the end of the primer (extract on to make 0-based for BED file)
                bed_line_hsp[2] = hsp.hit_end

                # Add the bed line for the hit to the bed file list
                bed_list.append(bed_line_hsp)

    # Construct lists to hold bed file names, and a dict to hold the number of hits for each primer pair
    blast_hit_beds = []
    exclusion_list = []
    primer_hits = {}

    # Split up the bed lines into each primer and write a bed file for each primer
    for primer_pair in primer_pairs.keys():
        # Find bed lines for primer pair
        primer_bed_lines = [line for line in bed_list
                            if primer_pairs[primer_pair][0] == line[3] or primer_pairs[primer_pair][1] == line[3]]

        # Record the number of hits in a genome by a set of primers
        primer_hits[primer_pair] = len(primer_bed_lines)

        # Write the bed files for each primer pair
        bed_name = f'{genome_name}~~{primer_pair}.bed'
        with open(bed_name, 'w') as bed_file:
            writer = csv.writer(bed_file, delimiter='\t')
            for line in primer_bed_lines:
                writer.writerow(line)
        bed_file.close()
        blast_hit_beds.append(bed_name)

        # Sort the bed file
        sorted_bed = bedtools.BedTool(bed_name).sort()
        sorted_bed.saveas(bed_name)

        # Check if primers are to be excluded, if then copy primer bed to be used later in bedtools intersect
        if not include_primers:
            copyfile(f'{genome_name}~~{primer_pair}.bed', f'{genome_name}~~{primer_pair}_primers.bed')
            exclusion_list.append(f'{genome_name}~~{primer_pair}_primers.bed')

    return blast_hit_beds, exclusion_list, primer_hits


def write_bed_from_list_of_primers(list_of_primers, bed_file_name):
    # Write out new .bed file with only the extended primers.
    # Join the items in the lists (primers) into a sting
    list_of_primers = [' '.join(primer) for primer in list_of_primers]

    # Join the primers into one large sting
    list_of_primers = '\n'.join(list_of_primers)

    # Construct BedTools object
    bed_of_primers = bedtools.BedTool(list_of_primers, from_string=True)

    # Sort the BedTools object
    bed_of_primers = bed_of_primers.sort()

    # Save the primers into a .bed file
    bed_of_primers.saveas(bed_file_name)


def primer_reach_contig_end_calc(genome_file, max_primer_dist, primer_contig_hits):
    # Construct dict to hold the indexing of the fasta file
    fai_dict = {}
    # Read and search the .fai file to get the edge of contigs
    with open(genome_file, 'r') as fai_file:
        for line in fai_file.readlines():
            line = line.split('\t')
            fai_dict[line[0]] = int(line[1])
    fai_file.close()

    # Extract all lines from dict holding Bed intervals for primer hits
    intervals = [interval for key in primer_contig_hits.keys() for interval in primer_contig_hits[key]]

    # Insert line identifier for interaction matrix
    [line.append(str(i)) for i, line in enumerate(intervals)]

    # Construct matrix to hold info on contig ends reached
    end_reached_matrix = [[0] * len(intervals) for _ in range(len(intervals))]

    # Go through each primer hit to see if it reached the end of the contig under the max distance
    for i, interval in enumerate(intervals):
        # Check if the 5' and 3' can be reached using the max distance
        if int(interval[1]) - max_primer_dist < 0:
            end_reached_matrix[i][0] = 1

        if int(interval[2]) + max_primer_dist > fai_dict[interval[0]]:
            end_reached_matrix[i][1] = 1

    # Find how many ends each primer can reach
    end_sums = [sum(ends) for ends in end_reached_matrix]

    # Count the number of primers that reach one or both ends
    end_reaches = [end_sums.count(2), end_sums.count(1)]

    return end_reaches, end_sums, end_reached_matrix, intervals


def examine_flanking_regions(primer_contig_hits, max_primer_dist, genome_file, bed_file_name=None):
    '''
    Function that takes a set of primers on contigs, and search for how they can be connected using a genome and length of contigs.
    :param primer_contig_hits: dict with each key being a contig and values being a list of seed sequence hit on contigs as lists
    :param max_primer_dist: int for the maximum distance allowed between seed sequenes
    :param genome_file: .fai (fasta-index) of a genome used as info on contig sizes
    :param bed_file_name: Used for saving updated BED file with coordinates
    :return: Evidence level for a set of seed seuqences that have hit the genome.
    '''
    # If no limits are given and only one contig is hit then no pairs can be predicted
    if max_primer_dist == 0:
        return 2

    # If limits are given then check overlap of primers, if they are only found on one contig
    elif len(primer_contig_hits) == 1 and max_primer_dist > 0:
        # Extract the intervals found for primers
        intervals = primer_contig_hits[list(primer_contig_hits.keys())[0]]

        # Construct an interaction matrix for the primers
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
            bed_interval = bed_interval.slop(b=max_primer_dist, g=genome_file)

            # Reinsert the newly slopped interval
            intervals_tmp[i] = list(bed_interval[0])

            # Construct BedTool object from all primer hits in intervals
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

            # Filter out overlapping primers in which the primers are the same
            interacting_ids = [line[5] for line in interacting_lines if len(set(line[4].split(','))) >= 2]
            interacting_primer_names = [line[4] for line in interacting_lines if len(set(line[4].split(','))) >= 2]

            # Insert interacting primers in matrix (only in upper part of matrix):
            if len(interacting_ids) > 0:
                for cluster in zip(interacting_ids, interacting_primer_names):
                    # split the pairs in the zipped tuple and find the combinations of names and ids
                    id_combinations = combinations(cluster[0].split(','), 2)
                    name_combinations = combinations(cluster[1].split(','), 2)

                    # Go through each pair of id and name combination to record valid hits
                    for pair in zip(id_combinations, name_combinations):
                        # Check that the combination of primers are not the same primer and that primer being examined is involved
                        if pair[1][0] != pair[1][1] and str(i) in pair[0]:
                            interaction_matrix[int(pair[0][0])][int(pair[0][1])] = 1


        # Summarise all rows of the interaction matrix
        num_interactions = [sum(row) for row in interaction_matrix]
        # Summarise col values to get a number of total interactions
        num_interactions = sum(num_interactions)

        if num_interactions == 0:
            return 1

        elif num_interactions == 1:
            # Check if primers reach end of contigs
            end_reaches, *_ = primer_reach_contig_end_calc(genome_file, max_primer_dist, primer_contig_hits)
            if sum(end_reaches) > 0:
                return 2

            interaction_indexes = []
            # Get index of the interaction between primers
            for i, interacts in enumerate(interaction_matrix):
                if 1 in interacts:
                    for j, interaction in enumerate(interaction_matrix[i]):
                        if interaction == 1:
                            interaction_indexes = [i, j]

            # Extract the primers that correspond to the identified indexes
            interaction_primers = [intervals[index][:4] for index in interaction_indexes]
            # Write the bed file from the list of primers
            write_bed_from_list_of_primers(interaction_primers, bed_file_name)
            return '5B'

        elif num_interactions > 1:
            return 2

    else:
        end_reaches, end_sums, end_reached_matrix, intervals = primer_reach_contig_end_calc(genome_file, max_primer_dist, primer_contig_hits)

        # See if no ends has been reached
        if end_reaches[0] == 0 and end_reaches[1] == 0:
            return 1

        # Check if one primer has reached one or both ends, while the other has reached no ends. # TODO - should this be handled differently? is it in line with evidence levels? Should the evidence level be changed for this?
        elif (end_reaches[0] == 1 and end_reaches[1] == 0) or (end_reaches[0] == 0 and end_reaches[1] == 1):
            return 6 # TODO Discuss the evidence level of this. - Maybe 2? # - test

        # Check if
        # two or more primers hit both ends of their contig (end_reaches[0] > 1) or
        # three or more primers all hit one end of their contig (end_reaches[1] > 2) or
        # one primer reaches two ends and two primers reaches one end each (end_reaches[0] == 1 and end_reaches[1] == 2)
        elif end_reaches[0] > 1 or end_reaches[1] > 2 or (end_reaches[0] == 1 and end_reaches[1] == 2):
            return 2

        # Check if one primer has reached both ends and one primer has reached one end
        elif end_reaches[1] == 1 and end_reaches[0] == 1:
            return 2

        # Check if two primers have reached one end.
        elif end_reaches[1] == 2:
            # Organise the bed and reprint with the usefull primers and adjusted coordinates to end
            primer_indexes = [index for index, reached in enumerate(end_sums) if reached == 1]

            # Initialise list to hold extended primers
            extended_primers = []

            for index in primer_indexes:
                # Get the primer extension hits
                primer_extension_hits = end_reached_matrix[index]

                # Get the primer in question, and remove the identifier at the end.
                extend_primer = intervals[index][:4]

                # Convert primer to a BedTools object
                extend_primer = bedtools.BedTool(' '.join(extend_primer), from_string=True)

                # Examine if primer should be extended in 5' end
                if primer_extension_hits[0] == 1:
                    extend_primer = extend_primer.slop(g=genome_file, l=max_primer_dist, r=0)

                # Examine if the primer should be extended in the 3' end
                elif primer_extension_hits[1] == 1:
                    extend_primer = extend_primer.slop(g=genome_file, r=max_primer_dist, l=0)

                else:
                    warnings.warn("A primer found to reach the end of a contig's "
                                  "end was not extended in any direction!\n"
                                  f"Please report this along with the genome: {genome_file},"
                                  f" the primer in question: {extend_primer}, "
                                  f"and the max distance set: {max_primer_dist}")
                    return "Something unexpected was found, please check warning log!"

                extended_primers.append(extend_primer[0])

            # Write out new .bed file with only the extended primers.
            write_bed_from_list_of_primers(extended_primers, bed_file_name)
            return '4B'

        else:
            # TODO - add what the tmp folder ends up being named.
            warnings.warn(f'Some unaccounted for constalation of primers hits '
                          f'was found to reach the end of contigs.\n'
                          f'Genome .fai file examined: {genome_file}.\n'
                          f'Please attach the .fai file (found in the tmp folder)'
                          f' when reporting this, along with the following: '
                          f'{intervals} and the max distance used: {max_primer_dist}')
            return "Something unexpected was found, please check warning log!"


def check_primers_placement(bed_files, primer_pairs, primer_hits, max_primer_dist, genome_file, file_type, tmp_folder):
    ''' Function to determine the placement of primers and how to search for connections between them
    Returns a genome file in the temporary directory and an evidence score for each bed file and its primer connections'''
    # Copy fasta to tmp folder to avoid .fai in input folder and possible clashed with permissions
    if file_type == 'fasta':
        tmp_genome = os.path.join(tmp_folder, genome_file.rsplit('/')[-1])
        copyfile(genome_file, tmp_genome)
        genome_file = tmp_genome

    # Produce genome index file using samtools
    samtools_faidx_cmd = SamtoolsFaidxCommandline(ref=genome_file)
    samtools_faidx_cmd()

    # Initialise the dict to hold the support for a primer hit
    primer_hit_support_dict = dict.fromkeys(primer_pairs)

    # Check each BED file
    for file in bed_files:
        # Strip file name to get name of genome
        primer_name = file.rsplit('.', 1)[0]
        primer_name = primer_name.rsplit('~~', 1)[-1]

        # Initialise dict to hold information on primer distribution over contigs
        primer_to_contig = {}

        # Check if more than one primers has hit the genome,
        # if, then check if they are same contig
        if primer_hits[primer_name] > 1:
            # Read Bed file
            bed_file = bedtools.BedTool(file)

            # Index lines in BED file and which contig seed sequences hit
            for line in bed_file:
                contig = line[0]
                try:
                    # Add the line to the contig
                    primer_to_contig[contig] += [list(line)]
                except KeyError:
                    # If not found add the contig and the line
                    primer_to_contig[contig] = [list(line)]

            # Find number of different primers that hit the contigs
            # Get all hits on contigs an unpack them into a single list of lists, then count the number of unique seed sequence hits
            primer_bed_lines = [primer_to_contig[key] for key in primer_to_contig.keys()]
            primer_bed_lines = [hit for contig in primer_bed_lines for hit in contig]
            uniq_primers = len(set([hit[3] for hit in primer_bed_lines]))

            # If more than one contig is hit and more than one seed sequence has hit then examine further,
            # else only one contig is hit
            if len(primer_to_contig.keys()) > 1 and uniq_primers > 1:
                # Examine if some of the primers are on the same contig:
                multi_hit_contigs = [contig for contig in primer_to_contig.keys() if len(primer_to_contig[contig]) >= 2]
                # Check if there are any contigs hit with multiple primers,
                # if, then examine them for overlaps between primers,
                # if not, examine across contigs
                if len(multi_hit_contigs) > 0:
                    # Examine all contigs that contain multiple primers for an overlap
                    for contig in multi_hit_contigs:
                        hit_contig = {contig: primer_to_contig[contig]}

                        return_value = examine_flanking_regions(hit_contig, max_primer_dist, f'{genome_file}.fai', file+'_5')

                        # Record the return value as evidence (3 = multiple overlaps)
                        # 2 or 4 means more examination is required
                        if return_value == 2:
                            primer_hit_support_dict[primer_name] = return_value
                        elif return_value == 1 or return_value == 4:
                            # Examine the remaining junctions to see if two primers can be found to connect across contigs
                            primer_hit_support_dict[primer_name] = examine_flanking_regions(primer_to_contig,
                                                                                            max_primer_dist,
                                                                                            f'{genome_file}.fai',
                                                                                            file)
                        elif return_value == '5B':
                            primer_hit_support_dict[primer_name] = return_value
                            alternative_return_value = examine_flanking_regions(primer_to_contig,
                                                                                max_primer_dist,
                                                                                f'{genome_file}.fai',
                                                                                file+'_6')

                            if alternative_return_value == '4B': # ADD alternative value 3?
                                primer_hit_support_dict[primer_name] = 2
                                os.remove(file+'_5')
                                try:
                                    os.remove(file + '_6')
                                except FileNotFoundError:
                                    pass

                            else:
                                primer_hit_support_dict[primer_name] = return_value
                                os.rename(file+'_5', file)

                        else:
                            sys.exit(100)
                else:
                    # Examine the primer hits across different contigs ***
                    return_value = examine_flanking_regions(primer_to_contig,
                                                            max_primer_dist,
                                                            f'{genome_file}.fai',
                                                            file)

                    if return_value == 1 and primer_hits[primer_name] == 2 and uniq_primers == 2:
                        primer_hit_support_dict[primer_name] = '4A'
                    else:
                        primer_hit_support_dict[primer_name] = return_value

            else:
                # Check if precisely two primers hit the contig and that the mates from the pair hit one time each
                if primer_hits[primer_name] == 2 and uniq_primers == 2:
                    # Score primer hit
                    primer_hit_support_dict[primer_name] = '5A' # TODO - POSSIBLY WRONG - WAS 7

                # Check if two unique primers hit, but one/both may have hit multiple times
                elif primer_hits[primer_name] > 2 and uniq_primers == 2:
                    # Examine primer hits and record returned score
                    primer_hit_support_dict[primer_name] = examine_flanking_regions(primer_to_contig,
                                                                                    max_primer_dist,
                                                                                    f'{genome_file}.fai')

                # Check that only one primer has hit, but it has hit multiple times
                elif uniq_primers == 1:
                    # Score primer hit
                    primer_hit_support_dict[primer_name] = 0
                else:
                    raise NotImplementedError(f'Some unaccounted for constalation of primers hits '
                                              f'was found to hit a single contig.\n'
                                              f'The primer pair in question is {primer_name} in genome {genome_file}.\n'
                                              f'Please report this along with the following: {primer_bed_lines}')
        else:
            primer_hit_support_dict[primer_name] = 0

    return genome_file, primer_hit_support_dict


def bed_merge_handling(blast_hit_beds, include_primers, exclude_primer_list, max_primer_dist, primer_evidence):
    # initialise list to hold merged bed file names:
    merged_bed_files = []

    # Process the bed file for each primer pair
    for i, bed_file in enumerate(blast_hit_beds):
        primer_hits = bedtools.BedTool(bed_file)

        # Get name of primer
        split_bed_name = bed_file.rsplit('/', 1)[-1]
        split_bed_name = split_bed_name.replace('.bed', '')
        primer_name = split_bed_name.rsplit('~~', 1)[-1]

        # Check if max distance between primers is set to unlimited (0)
        if max_primer_dist == 0:
            max_primer_dist = 9999999999

        # Merge the bed file and collapse the column containing the names and count the number of intervals collapsed
        if len(primer_hits) > 0:
            primer_hits = primer_hits.merge(c=[4, 4], o='collapse,count', d=max_primer_dist)
        else:
            primer_evidence[primer_name] = 0
            continue

        # Check if the evidence level for the primers should increase due to primers being merged.
        # Check that there is one line and that two primers have been merged and not just one primer hit
        if len(primer_hits) == 1 and int(primer_hits[0][4]) == 2:
            primer_evidence[primer_name] = '5B'

        # Evaluate if primers are to be excluded from the intervals
        # TODO - Can be problematic if the primer hits the end of a contig, then the primer will be removed. Either make the user aware that a primer is found on the edge of a contig and only extract sequences that have somthing between it and the contig break or other soluton.
        if not include_primers:
            # Load bed containing primers to be excluded
            exclusion_bed = bedtools.BedTool(exclude_primer_list[i])

            pre_deletion_intervals = len(primer_hits)

            # Remove the primer intervals
            primer_hits = primer_hits.subtract(exclusion_bed)

            if pre_deletion_intervals > len(primer_hits):
                primer_evidence[primer_name] = 4

        # Save the merged intervals
        # Construct name
        merged_bed_file = bed_file.rsplit('.', 1)[0]
        merged_bed_file = f'{merged_bed_file}_merged.bed'
        # Save
        # Check if there is any interval to report,
        # else notify by setting the return evidence level
        if len(primer_hits):
            primer_hits.saveas(merged_bed_file)
            # Add to list
            merged_bed_files.append(merged_bed_file)
        else:
            primer_evidence[primer_name] = 3

    return merged_bed_files, primer_evidence


def extract_seqs_n_annots(merged_bed_files, file_type, genome_file, annotation_file, tmp_folder, out_path, primer_pairs,
                          primer_evidence):
    # Initiate dict to hold number of annotations per interval
    annots_per_interval = dict.fromkeys(primer_pairs)
    # Initiate dict to hold primers that may neighbour a sequence break
    break_primers = {}
    # Initiate dict tot hold distances between primers
    inter_primer_dist = {}

    # Extract sequence from Fasta from all merged bed files
    for merged_bed in merged_bed_files:
        # Read Bed file
        merged_intervals = bedtools.BedTool(merged_bed)

        # Check if the bed file is empty indicate that no intervals are present by giving nan in annotations
        # There are no intervals when excluding primers and no primers are merged for specific primer set.
        if file_type == 'gff':
            primer_pair = merged_bed.rsplit('~~', 1)[-1]
            primer_pair = primer_pair.rsplit('_merged', 1)[0]
            if len(merged_intervals) == 0:
                annots_per_interval[primer_pair] = 'nan'
                continue

        # Go though each line of bed, really only necessary if there are breaks in genome.
        for interval in merged_intervals:
            # Identify primer pair for interval:
            primers = interval[3].split(',')
            # Check if primers are merged, or if single primer interval
            if len(primers) == 2:
                primer_pair_name = [primer_name for primer_name in primer_pairs if
                                    primers[0] in primer_pairs[primer_name] and primers[1] in primer_pairs[primer_name]]

                primer_pair_name = primer_pair_name[0]
            else:
                primer_pair_name = f'{primers[0]}_break'
                break_primers[primer_pair_name] = [primers[0], 'break']

            # Construct output name for extracted sequence
            genome_name = genome_file.rsplit('/', 1)[-1]
            genome_name = genome_name.rsplit('.', 1)[0]
            genome_name = genome_name.rsplit('_tmp', 1)[0]
            out_file_name = f'{genome_name}--{primer_pair_name}'
            output_genome = os.path.join(out_path, f'{out_file_name}.fasta')

            # Construct good fasta header
            fasta_header = f'{genome_name}_{primer_pair_name}'

            # Convert interval to string to be made into a BedTool object
            interval = f'{interval[0]} {int(interval[1])} {int(interval[2])} {fasta_header}'

            # Create BedTool object from line
            interval = bedtools.BedTool(interval, from_string=True)

            # Record the length of the segments between primers
            inter_primer_dist[primer_pair_name] = len(interval[0])

            # Extract the fasta sequence and save it in output folder
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
                annots_per_interval[primer_pair_name] = len(annot)

                # Start writing the output file, if an annotation if found
                if len(annot) > 0:
                    # Check if seed sequence reach a contig break
                    # If then get the primer pair to adjust the evidence level
                    if '_break' in primer_pair_name:
                        primer_pair_key = primer_pair_name
                        primer_pair_key = primer_pair_key.rsplit('_', 1)[0]
                        primer_pair_key = [key for key in primer_pairs.keys() if primer_pair_key in primer_pairs[key]][0]
                    else:
                        primer_pair_key = primer_pair_name

                    # Increase the evidence level for the seed sequence pair.
                    # Only if the evidence level is 4B or 5B and primer pair is valid (no breaks)
                    if primer_pair_key in primer_evidence.keys():
                        if primer_evidence[primer_pair_key] == '4B':
                            primer_evidence[primer_pair_key] = '4C'
                        elif primer_evidence[primer_pair_key] == '5B':
                            primer_evidence[primer_pair_key] = '5C'

                    # Construct name for output file
                    output_gff = os.path.join(out_path, f'{out_file_name}.gff')

                    # Find start coordinate for current merged interval
                    start_coordinate = int(interval[0][1])+1

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
        os.remove(f'{genome_file}.fai')
    else:
        os.remove(genome_file)
        os.remove(f'{genome_file}.fai')
        os.remove(annotation_file)

    return annots_per_interval, break_primers, primer_evidence, inter_primer_dist


def screen_genome_for_primers(genome_file, primer_pairs, primer_path, tmp_folder,
                              include_primers, file_type, annotation_file, out_path, max_primer_dist):
    # Clean the genome name for path, .gff and possible _tmp if gff is given
    # genome_name = genome_file.rsplit('/', 1)[1]
    genome_name = os.path.basename(genome_file)
    genome_name = genome_name.rsplit('.', 1)[0]
    genome_name = genome_name.rsplit('_tmp', 1)[0]

    # Concatenate the genome name and the path to the temporary folder
    genome_name = os.path.join(tmp_folder, genome_name)

    # Run blast with genome and insertion site sequences
    blast_xml_output = blast_insertion_site(primer_path, genome_file, f'{genome_name}_blast')

    # Construct a bedfile from blast output
    blast_hit_beds, exclude_primer_list, primer_hits = blast_out_to_sorted_bed(blast_xml_output,
                                                                               include_primers,
                                                                               genome_name,
                                                                               primer_pairs)

    # Examine the primer hits and try to find solutions when multiple primers hit at once
    genome_file, primer_evidence = check_primers_placement(blast_hit_beds, primer_pairs, primer_hits, max_primer_dist,
                                                           genome_file, file_type, tmp_folder)

    # sort and merge the bed files
    merged_bed_files, primer_evidence = bed_merge_handling(blast_hit_beds, include_primers, exclude_primer_list,
                                                           max_primer_dist, primer_evidence)

    # Extract sequences and annotations using merged intervals.
    annots_per_interval, break_primers, primer_evidence, inter_primer_dist = extract_seqs_n_annots(merged_bed_files,
                                                                                                   file_type,
                                                                                                   genome_file,
                                                                                                   annotation_file,
                                                                                                   tmp_folder,
                                                                                                   out_path,
                                                                                                   primer_pairs,
                                                                                                   primer_evidence)

    # clean up by removing blast xml output, and merged and primer hit beds from tmp folder:
    os.remove(blast_xml_output)
    [os.remove(file) for file in blast_hit_beds]
    [os.remove(file) for file in merged_bed_files]
    [os.remove(file) for file in exclude_primer_list]

    return primer_hits, annots_per_interval, genome_name, primer_evidence, break_primers, inter_primer_dist

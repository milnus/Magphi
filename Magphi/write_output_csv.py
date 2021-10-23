import csv
import os

# TODO - test?
def write_primer_hit_matrix(master_primer_hits, primer_pairs, out_path):
    with open(os.path.join(out_path, 'contig_hit_matrix.csv'), 'w') as out_file:
        header = ['genome'] + list(primer_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        for line in master_primer_hits.keys():
            writer.writerow(master_primer_hits[line])

        out_file.close()


def write_annotation_num_matrix(master_annotation_hits, primer_pairs, out_path):
    with open(os.path.join(out_path, 'annotation_num_matrix.csv'), 'w') as out_file:
        header = ['genome'] + list(primer_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        for line in master_annotation_hits.keys():
            writer.writerow(master_annotation_hits[line])

        out_file.close()


def write_primer_hit_evidence(master_primer_evidence, primer_pairs, out_path):
    with open(os.path.join(out_path, 'master_primer_evidence.csv'), 'w') as out_file:
        header = ['genome'] + list(primer_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        for line in master_primer_evidence.keys():
            writer.writerow(master_primer_evidence[line])

        out_file.close()


def write_inter_primer_dist(master_inter_primer_dist, primer_pairs, out_path):
    with open(os.path.join(out_path, 'inter_primer_distance.csv'), 'w') as out_file:
        header = ['genome'] + list(primer_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        for line in master_inter_primer_dist.keys():
            writer.writerow(master_inter_primer_dist[line])

        out_file.close()

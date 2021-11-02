import csv
import os


def write_seed_hit_matrix(master_seed_hits, seed_pairs, out_path):
    with open(os.path.join(out_path, 'contig_hit_matrix.csv'), 'w') as out_file:
        header = ['genome'] + list(seed_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        keys = list(master_seed_hits.keys())
        keys.sort()
        for line in keys:
            writer.writerow(master_seed_hits[line])

        out_file.close()


def write_annotation_num_matrix(master_annotation_hits, seed_pairs, out_path):
    with open(os.path.join(out_path, 'annotation_num_matrix.csv'), 'w') as out_file:
        header = ['genome'] + list(seed_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        keys = list(master_annotation_hits.keys())
        keys.sort()
        for line in keys:
            writer.writerow(master_annotation_hits[line])

        out_file.close()


def write_seed_hit_evidence(master_seed_evidence, seed_pairs, out_path):
    with open(os.path.join(out_path, 'master_seed_evidence.csv'), 'w') as out_file:
        header = ['genome'] + list(seed_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        keys = list(master_seed_evidence.keys())
        keys.sort()
        for line in keys:
            writer.writerow(master_seed_evidence[line])

        out_file.close()


def write_inter_seed_dist(master_inter_seed_dist, seed_pairs, out_path):
    with open(os.path.join(out_path, 'inter_seed_distance.csv'), 'w') as out_file:
        header = ['genome'] + list(seed_pairs.keys())
        writer = csv.DictWriter(out_file, fieldnames=header)

        # Write the field names or header row in file
        writer.writeheader()

        # Write the remaining lines
        keys = list(master_inter_seed_dist.keys())
        keys.sort()
        for line in keys:
            writer.writerow(master_inter_seed_dist[line])

        out_file.close()

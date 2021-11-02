import os
import warnings


def partition_outputs(seed_pairs, out_path, file_logger):
    """
    Function to distribute output files into appropriate output folders
    :param seed_pairs: Name of seed sequence pairs in a list
    :param out_path: Path to the output folder
    :param file_logger: Logger that outputs files to log
    :return:
    """
    file_logger.debug('Partitioning output files into seed folders')
    # Go through each seed sequence pair, construct sub output folder and populate with output files.
    for seed in seed_pairs:
        # List all output files in output folder
        out_folder_content = os.listdir(out_path)
        out_folder_content_paths = [os.path.join(out_path, file) for file in out_folder_content]

        # Construct folder for current seed to contain fasta and gff files for that seed
        seed_out_path = os.path.join(out_path, seed)
        try:
            os.mkdir(seed_out_path)
        except FileExistsError:
            file_logger.warning(f"Output folder for seed sequence pair: '{seed}' already exists")

        # find and move files into folder
        for i, file in enumerate(out_folder_content_paths):

            # Check that the marker of a output file is present (--)
            if '--' in file and ('.gff' in file or '.fasta' in file):
                # Isolate the seed from the file name
                seed_end = file.rsplit('--', 1)[-1]
                file_seed = seed_end.rsplit('.', 1)[0]

                # Check if break is present, if then remove the break and then narrow down the seed the file belongs to.
                if 'break' in file_seed:
                    file_seed = file_seed.rsplit('_', 1)[0]
                    seed_match = [file_seed != seed_name for seed_name in seed_pairs]

                    round_count = 1
                    while all(seed_match):
                        file_seed = file_seed[:-1]
                        seed_match = [file_seed != seed_name for seed_name in seed_pairs]

                        round_count += 1

                        if round_count == 100:
                            warnings.warn(f"Problem in placing outfile: '{file}' in an appropriate output folder for a seed.\n"
                                          f"Please report this!")
                            file_logger.warning(f"Problem in placing outfile: '{file}' in an appropriate output folder for a seed.\n"
                                                f"Please report this!")
                            continue

                # Check if a file in the output folder contain the seed name and is a fasta or gff file.
                # if seed == file_seed and (f'{seed}_' in seed_end or f'{seed}.' in seed_end):
                if seed == file_seed:
                    # Replace the double hyphen with a single.
                    file_new = out_folder_content[i].replace('--', '-')
                    # Move the file into the seed folder
                    os.rename(file, os.path.join(seed_out_path, file_new))


def write_paired_seeds(seed_pairs, out_path):
    with open(os.path.join(out_path, 'seed_pairing.tsv'), 'w') as seed_file:
        for seed_key in seed_pairs.keys():
            # Fetch the seeds
            seeds = seed_pairs[seed_key]
            # Write the line
            seed_file.write(f'{seed_key}\t{seeds[0]}\t{seeds[1]}\n')
        seed_file.close()

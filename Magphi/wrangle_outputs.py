import os
import warnings


def partition_outputs(primer_pairs, out_path):
    # Go through each primer set, construct folder and populate the folder.
    for primer in primer_pairs:
        # List all output files in output folder
        out_folder_content = os.listdir(out_path)
        out_folder_content_paths = [os.path.join(out_path, file) for file in out_folder_content]

        # Construct folder for current primer to contain fasta and gff files for that primer
        primer_out_path = os.path.join(out_path, primer)
        try:
            os.mkdir(primer_out_path)
        except FileExistsError:
            Warning(f"Output folder for primer pair: {primer} already exists")

        # find and move files into folder
        for i, file in enumerate(out_folder_content_paths):

            # Check that the marker of a output file is present (--)
            if '--' in file and ('.gff' in file or '.fasta' in file):
                # Isolate the primer from the file name
                primer_end = file.rsplit('--', 1)[-1]
                file_primer = primer_end.rsplit('.', 1)[0]

                # Check if break is present, if then remove the break and then narrow down the primer the file belongs to.
                if 'break' in file_primer:
                    file_primer = file_primer.rsplit('_', 1)[0]
                    primer_match = [file_primer != primer_name for primer_name in primer_pairs]

                    round_count = 1
                    while all(primer_match):
                        file_primer = file_primer[:-1]
                        primer_match = [file_primer != primer_name for primer_name in primer_pairs]

                        round_count += 1

                        if round_count == 100:
                            warnings.warn(f'Problem in placing outfile: {file} in an appropriate output folder for a primer.\n'
                                          f'Please report this!')
                            continue

                # Check if a file in the output folder contain the primer name and is a fasta or gff file.
                # if primer == file_primer and (f'{primer}_' in primer_end or f'{primer}.' in primer_end):
                if primer == file_primer:
                    # Replace the double hyphen with a single.
                    file_new = out_folder_content[i].replace('--', '-')
                    # Move the file into the primer folder
                    os.rename(file, os.path.join(primer_out_path, file_new))


def write_paired_primers(primer_pairs, out_path):
    with open(os.path.join(out_path, 'primer_pairing.tsv'), 'w') as primer_file:
        for primer_key in primer_pairs.keys():
            # Fetch the primers
            primers = primer_pairs[primer_key]
            # Write the line
            primer_file.write(f'{primer_key}\t{primers[0]}\t{primers[1]}\n')
        primer_file.close()

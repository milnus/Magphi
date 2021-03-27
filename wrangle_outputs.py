import os


def partition_outputs(primer_pairs, out_path):
    # Go through each primer that set, construct folder and populate the folder.
    for primer in primer_pairs:
        # List all output files in output folder
        out_folder_content = os.listdir(out_path)
        out_folder_content_paths = [os.path.join(out_path, file) for file in out_folder_content]

        # Construct folder for each primer to contain fasta and gff files for that primer
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
                file_primer = file.rsplit('.', 1)[0]
                file_primer = file_primer.rsplit('--', 1)[-1]
                # Check if break is present
                if 'break' in file_primer:
                    file_primer = file_primer.rsplit('_', 2)[0]

                # Check if a file in the output folder contain the primer name and is a fasta or gff file.
                if primer in file_primer:
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

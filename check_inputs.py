def check_if_fasta(input_files):
    # Construct list to hold if files are fasta
    is_inputs_fasta = [False]*len(input_files)

    # Go through each file and check if it is a fasta genome
    for i, file in enumerate(input_files):
        with open(file, 'r') as in_file:
            if '>' in in_file.readline():
                is_inputs_fasta[i] = True

    if all(is_inputs_fasta):
        # TODO - make verbose controlled
        print(f'all files were found to be fasta files')
        return 'fasta'

    elif any(is_inputs_fasta):
        print('UPS! Something went wrong!\n'
              'Input files found to be of mixed type! '
              'Please check that all your input files are either Fasta or GFF, not mixed')
        exit(code=1)

    else:
        return None


def check_if_gff(input_files):
    is_inputs_gff = [False] * len(input_files)
    for i, file in enumerate(input_files):
        with open(file, 'r') as in_file:
            if '##gff-version 3' in in_file.readline():
                for line in in_file.readlines():
                    if '##FASTA' in line:
                        is_inputs_gff[i] = True
                # Check that a genome has been found in file
                if is_inputs_gff[i] is False:
                    print(f'UPS! {file} does seem to be a GFF3 version, but not one that contain the genome '
                          f'following a ##FASTA line - Please check the file')
                    exit(code=1)

    if all(is_inputs_gff):
        # TODO - make verbose controlled
        print(f'all files were found to be GFF3 files')
        return 'gff'
    else:
        return None


def check_inputs(input_files):
    file_type = check_if_fasta(input_files)

    # Check if input files are GFF3
    if file_type is None:
        file_type = check_if_gff(input_files)

    # If file_type is still none, exit as input files are not recognized
    if file_type is None:
        print('UPS! Something went wrong!\n'
              'The given input files could not be recognised as either Fasta or GFF3 with a genome.')
        exit(code=1)

    return file_type

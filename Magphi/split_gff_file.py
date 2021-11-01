import os
import concurrent.futures
import gzip


def split_single_gff(gff, tmp_folder, is_input_gzipped):
    # Trim the name/path of the input file
    base_name = gff.split('/')[-1]
    base_name = base_name.rsplit('.gz', 1)[0]
    base_name = base_name.rsplit('.', 1)[0]

    # Add the path to the temporary folder
    outfolder_basename = os.path.join(tmp_folder, base_name)

    # Construct the names of the temporary files
    tmp_fasta_name = f'{outfolder_basename}_tmp.fasta'
    tmp_gff_name = f'{outfolder_basename}_tmp.gff'

    # Initialise the temporary files for the genome and annotations
    tmp_gff = open(tmp_gff_name, 'w')
    tmp_fasta = open(tmp_fasta_name, 'w')

    # Open the input gff file
    if is_input_gzipped:
        gff_file = gzip.open(gff, 'rt')
    else:
        gff_file = open(gff, 'r')
    # Indicate whether fasta file is found
    fasta_found = False

    # Go through all line of gff file
    for line in gff_file.readlines():
        # Check if genome has been reach if, then skip line
        if '##FASTA' in line:
            fasta_found = True
            continue

        # Check if genome has been reach if, then record line in genome,
        # else then record the as annotation line
        if fasta_found:
            tmp_fasta.write(line)
        else:
            tmp_gff.write(line)

    # Close files before opening the next
    gff_file.close()
    tmp_gff.close()
    tmp_fasta.close()

    return tmp_fasta_name, tmp_gff_name


def split_gff_files(input_files, tmp_folder, is_input_gzipped):
    """ Function to go through each input GFF3 file and split it into temporary files containing the gff annotations
    and the genome, to make them compatible with BLAST and bedtools"""

    # initialise lists to hold the name of the fasta and annotations files
    genome_files = []
    annotation_files = []

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = [executor.submit(split_single_gff, gff, tmp_folder, is_input_gzipped) for gff in input_files]

        for f in concurrent.futures.as_completed(results):
            tmp_fasta, tmp_gff = f.result()

            genome_files.append(tmp_fasta)
            annotation_files.append(tmp_gff)

    return genome_files, annotation_files

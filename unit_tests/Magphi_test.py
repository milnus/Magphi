'''
Unit tests for Magphi.

Usage: python -m unittest -v Magphi_test
'''

import unittest
import os
import json
from shutil import copyfile
import logging

from Magphi import check_inputs
from Magphi import split_gff_file
from Magphi import seed_handling
from Magphi import search_insertion_sites
from Magphi import wrangle_outputs
from Magphi import write_output_csv
from Magphi import exit_with_error
# pylint: disable=E1133

from io import StringIO
# pylint: disable=no-name-in-module

# Move to folder with mock input files. First try Github structure, then try pulled repository structure
try:
    os.chdir('/Magphi/unit_tests/unit_test_data/')
except FileNotFoundError:
    os.chdir('../unit_test_data/')


class TestExitWithError(unittest.TestCase):
    def test_exit_w_tmp_folder_deletion(self):
        ''' Test the exit function is able to remove the temporary folder '''
        tmp_folder = 'TestExitWithError/tmp_folder'
        tmp_folder_copy = 'TestExitWithError/tmp_folder_copy'
        os.mkdir(tmp_folder_copy)

        tmp_files = os.listdir(tmp_folder)
        for file in tmp_files:
            copyfile(os.path.join(tmp_folder, file), os.path.join(tmp_folder_copy, file))

        with self.assertRaises(SystemExit):
            exit_with_error.exit_with_error(exit_status=2, message='test msg', tmp_folder=tmp_folder)

        os.rename(tmp_folder_copy, tmp_folder)


class TestFileRecognition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_fasta_recognition(self):
        ''' test the recognition of fasta files '''
        path = 'TestFileRecognition/Fasta_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files, self.logger, False)

        self.assertEqual('fasta', file_type)

    def test_gzipped_fasta_recognition(self):
        ''' test the recognition of fasta files '''
        path = 'TestFileRecognition/All_gzipped_fasta'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files, self.logger, True)

        self.assertEqual('fasta', file_type)

    def test_none_fasta_recognition(self):
        ''' test that gff files are not recognised as fasta files '''
        path = 'TestFileRecognition/Gff3_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files, self.logger, False)

        self.assertEqual(None, file_type)

    def test_file_with_single_greater_than(self):
        ''' test that input files woth just a single '>' are not recognised as valid '''
        files = ['TestFileRecognition/single_greater_than.fasta']

        with self.assertRaises(SystemExit):
            check_inputs.check_if_fasta(files, self.logger, False)

    def test_file_with_new_line_in_sequence(self):
        ''' test that input files woth just a single '>' are not recognised as valid '''
        files = ['TestFileRecognition/new_line_in_seqeunce.fasta']

        with self.assertRaises(SystemExit):
            check_inputs.check_if_fasta(files, self.logger, False)

    def test_mixed_gff_and_fasta_recognition(self):
        ''' test that a mix of fasta and gff files results in exiting Magphi with an error '''
        path = 'TestFileRecognition/Mixed_gff_and_fasta'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_fasta(files, self.logger, False)

    def test_fasta_and_random_text_recognition(self):
        ''' test that a mix of fasta and random text files results in exiting Magphi with an error '''
        path = 'TestFileRecognition/Mixed_fasta_and_text'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_fasta(files, self.logger, False)

    def test_gff_recognition(self):
        ''' test that gff files with an attached genome are recognised correctly '''
        path = 'TestFileRecognition/Gff3_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_gff(files, self.logger, False)

        self.assertEqual('gff', file_type)

    def test_gzipped_gff_recognition(self):
        ''' test that gff files with an attached genome are recognised correctly '''
        path = 'TestFileRecognition/All_gzipped_gff'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_gff(files, self.logger, True)

        self.assertEqual('gff', file_type)

    def test_none_gff_recognition(self):
        ''' test that fasta files are not recognised as gff '''
        path = 'TestFileRecognition/Fasta_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_gff(files, self.logger, False)

        self.assertEqual(None, file_type)

    def test_gff_missing_genome_recognition(self):
        ''' test that gff files without a genomes attached exits with an error '''
        path = 'TestFileRecognition/Gff3_without_genome_attached'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]
        with self.assertRaises(SystemExit):
            check_inputs.check_if_gff(files, self.logger, False)

    def test_gff_and_random_text_recognition(self):
        ''' test that a mix of GFF3 and random text files results in exiting Magphi with an error '''
        path = 'TestFileRecognition/Mixed_gff_and_text'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_gff(files, self.logger, False)

    def test_not_incompatible_recognition(self):
        ''' test that a text file not being a Fasta or GFF3 files results in an error '''
        files = ['TestFileRecognition/Mixed_miscellaneous_files/Random_text.txt']

        with self.assertRaises(SystemExit):
            check_inputs.check_inputs(files, self.logger)

    def test_empty_file(self):
        ''' Test that an empty file results in an error '''
        files = ['TestFileRecognition/Mixed_miscellaneous_files/empty_file.txt']

        with self.assertRaises(SystemExit):
            check_inputs.check_inputs(files, self.logger)

    def test_gzipped_files(self):
        ''' Test that input of gzipped files are recognised as such '''
        files = ['TestFileRecognition/All_gzipped_fasta/GCA_005163865.fna.gz', 'TestFileRecognition/All_gzipped_fasta/GCA_900475985.fna.gz']

        self.assertEqual(True, check_inputs.check_if_gzip(files, self.logger))

    def test_mixed_gzipped_and_none_compressed_files(self):
        ''' Test that input of gzipped files are recognised as such '''
        files = ['TestFileRecognition/Mixed_gzipped/GCA_005163865.fna.gz', 'TestFileRecognition/Mixed_gzipped/GCA_900475985.fna']

        with self.assertRaises(SystemExit):
            check_inputs.check_if_gzip(files, self.logger)

    def test_non_gzipped_files(self):
        ''' Test that input of gzipped files are recognised as such '''
        files = ['TestFileRecognition/Fasta_files/GCA_005163865.fna']

        self.assertEqual(False, check_inputs.check_if_gzip(files, self.logger))


class TestSplittingGff(unittest.TestCase):
    def test_gff_split_single_file(self):
        ''' test the function that splits a gff file into annotations and genome. Assess the number of lines in output '''
        path = os.getcwd()
        file = os.path.join(path, 'TestSplittingGff/minimized.gff')

        # Split the test file
        genome, annotation = split_gff_file.split_single_gff(file, path, is_input_gzipped=False)

        # read the now divided genome and annotations and get the number of lines
        open_genome = open(genome, 'r')
        open_annotation = open(annotation, 'r')
        genome_file_length = len(open_genome.readlines())
        annotation_file_length = len(open_annotation.readlines())

        # Close files again
        open_genome.close()
        open_annotation.close()

        # Test if the files contain the number of expected lines.
        self.assertEqual(10, genome_file_length)
        self.assertEqual(5, annotation_file_length)

        # remove the genome and annotations.
        os.remove(genome)
        os.remove(annotation)

    def test_gff_split_single_gzipped_file(self):
        ''' test the function that splits a gff file into annotations and genome. Assess the number of lines in output '''
        path = os.getcwd()
        file = os.path.join(path, 'TestSplittingGff/minimized_gzipped.gff.gz')

        # Split the test file
        genome, annotation = split_gff_file.split_single_gff(file, path, is_input_gzipped=True)

        # read the now divided genome and annotations and get the number of lines
        open_genome = open(genome, 'r')
        open_annotation = open(annotation, 'r')
        genome_file_length = len(open_genome.readlines())
        annotation_file_length = len(open_annotation.readlines())

        # Close files again
        open_genome.close()
        open_annotation.close()

        # Test if the files contain the number of expected lines.
        self.assertEqual(10, genome_file_length)
        self.assertEqual(5, annotation_file_length)

        # remove the genome and annotations.
        os.remove(genome)
        os.remove(annotation)


class TestseedFunctions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_uneven_seed_number(self):
        ''' test that program exits if an uneven number of seeds is given, as this can not be made into a number of sets '''
        with self.assertRaises(SystemExit):
            seed_handling.check_number_n_names_of_seeds('TestPrimerFunctions/Uneven_number_primers.txt', self.logger)

    def test_identical_seed_names(self):
        ''' Test that giving seeds with the exact same name will result in an exit of the program '''
        with self.assertRaises(SystemExit):
            seed_handling.check_number_n_names_of_seeds('TestPrimerFunctions/Same_name_primers.txt', self.logger)

    def test_extracting_seed_info(self):
        ''' Test that giving seeds with the exact same name will result in an exit of the program '''
        num, names = seed_handling.check_number_n_names_of_seeds('TestPrimerFunctions/Good_seed_file.txt', self.logger)

        self.assertEqual(4, num)
        self.assertEqual(['primer_1_1', 'primer_1_2', 'primer_2_1', 'primer_2_2'], names)

    def test_correct_seed_pairing(self):
        ''' test that a file with correctly named seeds can be paired as expected '''
        seed_names = ['D_1', 'D_2', 'mutsD_1', 'mutsD_2']
        seed_pairs = seed_handling.construct_pair_seeds(seed_names, self.logger)

        expected_names = {'D': ['D_1', 'D_2'], 'mutsD': ['mutsD_1', 'mutsD_2']}

        # Sort dicts to make them same order
        seed_pairs = [seed_pairs[key].sort() for key in seed_pairs]
        expected_names = [expected_names[key].sort() for key in expected_names]

        self.assertEqual(expected_names , seed_pairs)

    def test_non_matching_seed_names(self):
        ''' Test that giving seeds that can not be matched by name makes the program exit '''
        seed_names = ['D_2', 'B_1', 'Z_1', 'A_1']
        with self.assertRaises(SystemExit):
            seed_handling.construct_pair_seeds(seed_names, self.logger)



# TODO - test blast function?


# TODO - test blast_out_to_sorted_bed function - use an input file of blast xml output - use two sets of seeds
#  Test both inclusion and exclution of seeds. - Andrew's responsitibily.
#   - We want to test that a blast output is converted correctly to Bed format.
#   1. Produce mock fasta to blast against. (Should have known sites that seeds match. Maybe repeat single Base or gap with seeds being unique)
#   2. Produce mock seeds
#   3. Blast mock seeds against mock fasta using similar settings as Magphi
#   4. Manually curate the positions are as expected
#   5. Manually determine the expected bed file information
#   6. Convert expected bed file format into .json or staight python code to be used for assertion.
#   7. write test.
#   8. Run


class TestseedsPlacement(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Construct mock logger
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def test_single_seed_single_hit_low_max_dist(self):
        ''' Test that a single seed sequence hit returns the correct evidence level '''
        bed_files = ['TestPrimersPlacement/single_contig_1200N~~single_primer.bed']
        seed_pairs = {'single_primer': ['single_primer_1', 'single_primer_2']}
        seed_hits = {'single_primer': 1}
        max_seed_dist = 1
        genome_file = 'TestPrimersPlacement/single_contig_1200N.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/single_contig_1200N.fasta.fai')

        evidence_level_return = flanking_return[1]['single_primer']
        self.assertEqual(0, evidence_level_return)

    def test_single_seed_single_hit_large_max_dist(self):
        ''' Test that a single seed sequence hit returns the correct evidence level '''
        bed_files = ['TestPrimersPlacement/single_contig_1200N~~single_primer.bed']
        seed_pairs = {'single_primer': ['single_primer_1', 'single_primer_2']}
        seed_hits = {'single_primer': 1}
        max_seed_dist = 10000
        genome_file = 'TestPrimersPlacement/single_contig_1200N.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/single_contig_1200N.fasta.fai')

        evidence_level_return = flanking_return[1]['single_primer']
        self.assertEqual(0, evidence_level_return)

    def test_single_seed_multiple_hit_same_contig_no_overlap(self):
        ''' Test that a single seed of a pair, hitting a single contig multiple times results in a correct evidence level '''
        bed_files = ['TestPrimersPlacement/single_contig_1200N~~primer.bed']
        seed_pairs = {'primer': ['primer_1', 'primer_2']}
        seed_hits = {'primer': 2}
        max_seed_dist = 1
        genome_file = 'TestPrimersPlacement/single_contig_1200N.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/single_contig_1200N.fasta.fai')

        evidence_level_return = flanking_return[1]['primer']
        self.assertEqual(0, evidence_level_return)

    def test_single_seed_multiple_hit_same_contig_w_overlap(self):
        ''' Test that a single seed hitting a single contig multiple times with large max distance gives the correct evidence level '''
        bed_files = ['TestPrimersPlacement/single_contig_1200N~~primer.bed']
        seed_pairs = {'primer': ['primer_1', 'primer_2']}
        seed_hits = {'primer': 2}
        max_seed_dist = 500
        genome_file = 'TestPrimersPlacement/single_contig_1200N.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/single_contig_1200N.fasta.fai')

        evidence_level_return = flanking_return[1]['primer']
        self.assertEqual(0, evidence_level_return)

    def test_single_seed_multiple_hit_multiple_contigs_no_overlap(self):
        ''' Test that a single seed hitting multiple contigs with no overlap results in the right evidence level '''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_same.bed']
        seed_pairs = {'primer_same': ['primer_same_1', 'primer_same_2']}
        seed_hits = {'primer_same': 2}
        max_seed_dist = 1
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        evidence_level_return = flanking_return[1]['primer_same']
        self.assertEqual(0, evidence_level_return)

    def test_single_seed_multiple_hit_multiple_contigs_with_overlap(self):
        ''' Test that a single seed hitting multiple contigs with overlap results in the right evidence level '''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_same.bed']
        seed_pairs = {'primer_same': ['primer_same_1', 'primer_same_2']}
        seed_hits = {'primer_same': 2}
        max_seed_dist = 1000
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        evidence_level_return = flanking_return[1]['primer_same']
        self.assertEqual(0, evidence_level_return)

    def test_multiple_seeds_multiple_hit_single_contig(self):
        ''' Test that both seed sequences hit once on one contig with no overlap and that the evidence level is correct '''
        bed_files = ['TestPrimersPlacement/single_contig_1200N~~primer_different.bed']
        seed_pairs = {'primer_different': ['primer_different_1', 'primer_different_2']}
        seed_hits = {'primer_different': 2}
        max_seed_dist = 1
        genome_file = 'TestPrimersPlacement/single_contig_1200N.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/single_contig_1200N.fasta.fai')

        evidence_level_return = flanking_return[1]['primer_different']
        self.assertEqual('5A', evidence_level_return)

    def test_multiple_hits_two_seeds_single_contig(self):
        ''' Test that two unique seed sequence hitting the same contig multiple times returns the correct evidence level '''
        bed_files = ['TestPrimersPlacement/single_contig_1200N~~primer_multi_different.bed']
        seed_pairs = {'primer_multi_different': ['primer_multi_different_1', 'primer_multi_different_2']}
        seed_hits = {'primer_multi_different': 3}
        max_seed_dist = 1
        genome_file = 'TestPrimersPlacement/single_contig_1200N.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)
        # remove .fai file
        os.remove('TestPrimersPlacement/single_contig_1200N.fasta.fai')

        evidence_level_return = flanking_return[1]['primer_multi_different']
        self.assertEqual(1, evidence_level_return)

    def test_multiple_hits_multiple_contigs_inter_contig_connect(self):
        ''' Test the outcome with two seed sequences that can connect on same contig, but not across contigs
        test both the altered bed file returned and the evidence level'''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_close_placement.bed']
        seed_pairs = {'primer_close_placement': ['primer_close_placement_1', 'primer_close_placement_2']}
        seed_hits = {'primer_close_placement': 2}
        max_seed_dist = 51
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'


        # Copy input bed file as it gets altered
        copyfile(bed_files[0], bed_files[0] + 'original')

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)

        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        # Check altered input file
        with open(bed_files[0], 'r') as altered_file:
            self.assertEqual(['Contig_1\t100\t300\tprimer_close_placement_1\n', 'Contig_1\t350\t550\tprimer_close_placement_2\n'],
                             altered_file.readlines())

        # Copy back input file
        os.rename(bed_files[0] + 'original', bed_files[0])
        evidence_level_return = flanking_return[1]['primer_close_placement']
        self.assertEqual('5B', evidence_level_return) #***

    def test_multiple_hits_multiple_contigs_cross_contig_connect(self):
        ''' Test the outcome with two unique seed sequences that can connect on same contig and across contigs'''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_close_placement.bed']
        seed_pairs = {'primer_close_placement': ['primer_close_placement_1', 'primer_close_placement_2']}
        seed_hits = {'primer_close_placement': 2}
        max_seed_dist = 101
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'


        # Copy input bed file as it gets altered
        copyfile(bed_files[0], bed_files[0]+'original')

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)

        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        # Copy back input file
        evidence_level_return = flanking_return[1]['primer_close_placement']
        self.assertEqual(2, evidence_level_return)

    def test_multiple_hits_multiple_contigs_cross_contig_reach(self):
        ''' Test the outcome when two unique seed sequences can connect across contigs when inter contig connection is not allowed by max distance'''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_long_placement.bed']
        seed_pairs = {'primer_long_placement': ['primer_long_placement_1', 'primer_long_placement_2']}
        seed_hits = {'primer_long_placement': 2}
        max_seed_dist = 101
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'


        # Copy input bed file as it gets altered
        copyfile(bed_files[0], bed_files[0] + 'original')

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)

        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        # Check altered input file
        with open(bed_files[0], 'r') as altered_file:
            self.assertEqual(
                ['Contig_1\t0\t300\tprimer_long_placement_1\n', 'Contig_2\t0\t75\tprimer_long_placement_2\n'],
                altered_file.readlines())


        # Copy back input file
        os.rename(bed_files[0] + 'original', bed_files[0])
        evidence_level_return = flanking_return[1]['primer_long_placement']
        self.assertEqual('4B', evidence_level_return)

    def test_multiple_hits_multiple_contigs_multi_overlap_long(self):
        ''' Test the outcome when two unique seed sequences can connect on same contig and across contigs'''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_long_placement.bed']
        seed_pairs = {'primer_long_placement': ['primer_long_placement_1', 'primer_long_placement_2']}
        seed_hits = {'primer_long_placement': 2}
        max_seed_dist = 401
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'

        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)

        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        evidence_level_return = flanking_return[1]['primer_long_placement']
        self.assertEqual(2, evidence_level_return)

    def test_multiple_hits_multiple_contigs_end_overlap_short(self):
        ''' Test the outcome with two unique seed sequences can connect across the contig and across the ends of the same contig'''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_short_placement.bed']
        seed_pairs = {'primer_short_placement': ['primer_short_placement_1', 'primer_short_placement_2']}
        seed_hits = {'primer_short_placement': 2}
        max_seed_dist = 126
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'


        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)

        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        evidence_level_return = flanking_return[1]['primer_short_placement']
        self.assertEqual(2, evidence_level_return)

    def test_multiple_hits_multiple_contigs_multi_overlap_short(self):
        ''' Test the outcome with two unique seed sequences can connect on same contig and across different contigs'''
        bed_files = ['TestPrimersPlacement/double_contig~~primer_short_placement.bed']
        seed_pairs = {'primer_short_placement': ['primer_short_placement_1', 'primer_short_placement_2']}
        seed_hits = {'primer_short_placement': 2}
        max_seed_dist = 1000
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'


        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)

        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        evidence_level_return = flanking_return[1]['primer_short_placement']
        self.assertEqual(2, evidence_level_return)

    def test_sinlge_hits_multiple_contigs_no_overlap(self):
        ''' Test the outcome with only two unique seed sequences can connect across different contigs'''
        bed_files = ['TestPrimersPlacement/double_contig~~single_primers_across_contigs.bed']
        seed_pairs = {'single_primers_across_contigs': ['single_primers_across_contigs_1', 'single_primers_across_contigs_2']}
        seed_hits = {'single_primers_across_contigs': 2}
        max_seed_dist = 1
        genome_file = 'TestPrimersPlacement/double_contig.fasta'
        file_type = 'fasta'
        tmp_folder = 'TestPrimersPlacement'


        flanking_return = search_insertion_sites.check_seeds_placement(bed_files=bed_files,
                                                                         seed_pairs=seed_pairs,
                                                                         seed_hits=seed_hits,
                                                                         max_seed_dist=max_seed_dist,
                                                                         genome_file=genome_file,
                                                                         file_type=file_type,
                                                                         tmp_folder=tmp_folder)

        os.remove('TestPrimersPlacement/double_contig.fasta.fai')

        evidence_level_return = flanking_return[1]['single_primers_across_contigs']
        self.assertEqual('4A', evidence_level_return)


class TestseedReachContigEndCalculation(unittest.TestCase):
    def test_no_ends_reached(self):
        ''' Test for a situation where no ends of contigs are reached by seeds '''
        genome_file_fai = 'TestPrimerReachContigEndCalculation/single_contig_1200N.fasta.fai'
        max_distance = 1
        seed_contig_hits = {'Contig_1': [['Contig_1', 500, 600, 'Primer_1'], ['Contig_1', 600, 700, 'Primer_2']]}

        end_reaches, end_sums, end_reached_matrix, intervals = \
            search_insertion_sites.seed_reach_contig_end_calc(genome_file_fai, max_distance, seed_contig_hits)

        self.assertEqual([0, 0], end_reaches)
        self.assertEqual([0, 0], end_sums)
        self.assertEqual([[0, 0], [0, 0]], end_reached_matrix)
        self.assertEqual([['Contig_1', 500, 600, 'Primer_1', '0'], ['Contig_1', 600, 700, 'Primer_2', '1']], intervals)

    def test_single_3_prime_end_reached(self):
        ''' test for a single contig that reach only the 3' of a contig. '''
        genome_file_fai = 'TestPrimerReachContigEndCalculation/single_contig_1200N.fasta.fai'
        max_distance = 500
        seed_contig_hits = {'Contig_1': [['Contig_1', 500, 600, 'Primer_1'], ['Contig_1', 600, 700, 'Primer_2']]}

        end_reaches, end_sums, end_reached_matrix, intervals = \
            search_insertion_sites.seed_reach_contig_end_calc(genome_file_fai, max_distance, seed_contig_hits)

        self.assertEqual([0, 1], end_reaches)
        self.assertEqual([0, 1], end_sums)
        self.assertEqual([[0, 0], [0, 1]], end_reached_matrix)
        self.assertEqual([['Contig_1', 500, 600, 'Primer_1', '0'], ['Contig_1', 600, 700, 'Primer_2', '1']], intervals)

    def test_single_5_prime_end_reached(self):
        ''' Test for a single seed sequence that each the 5' of a contig. '''
        genome_file_fai = 'TestPrimerReachContigEndCalculation/single_contig_1200N.fasta.fai'
        max_distance = 451
        seed_contig_hits = {'Contig_1': [['Contig_1', 450, 600, 'Primer_1'], ['Contig_1', 600, 650, 'Primer_2']]}

        end_reaches, end_sums, end_reached_matrix, intervals = \
            search_insertion_sites.seed_reach_contig_end_calc(genome_file_fai, max_distance, seed_contig_hits)

        self.assertEqual([0, 1], end_reaches)
        self.assertEqual([1, 0], end_sums)
        self.assertEqual([[1, 0], [0, 0]], end_reached_matrix)
        self.assertEqual([['Contig_1', 450, 600, 'Primer_1', '0'], ['Contig_1', 600, 650, 'Primer_2', '1']], intervals)

    def test_single_seed_sequence_reach_both_ends(self):
        ''' Test for a seed sequence that reach both ends of a contig '''
        genome_file_fai = 'TestPrimerReachContigEndCalculation/single_contig_1200N.fasta.fai'
        max_distance = 451
        seed_contig_hits = {'Contig_1': [['Contig_1', 450, 750, 'Primer_1'], ['Contig_1', 600, 650, 'Primer_2']]}

        end_reaches, end_sums, end_reached_matrix, intervals = \
            search_insertion_sites.seed_reach_contig_end_calc(genome_file_fai, max_distance, seed_contig_hits)

        self.assertEqual([1, 0], end_reaches)
        self.assertEqual([2, 0], end_sums)
        self.assertEqual([[1, 1], [0, 0]], end_reached_matrix)
        self.assertEqual([['Contig_1', 450, 750, 'Primer_1', '0'], ['Contig_1', 600, 650, 'Primer_2', '1']], intervals)

    def test_single_seed_sequence_reach_both_ends_n_seed_seqeunce_with_one_reach(self):
        ''' Test case of seed sequences that reach differing number of ends of contigs. '''
        genome_file_fai = 'TestPrimerReachContigEndCalculation/single_contig_1200N.fasta.fai'
        max_distance = 451
        seed_contig_hits = {'Contig_1': [['Contig_1', 450, 750, 'Primer_1'], ['Contig_1', 600, 750, 'Primer_2']]}

        end_reaches, end_sums, end_reached_matrix, intervals = \
            search_insertion_sites.seed_reach_contig_end_calc(genome_file_fai, max_distance, seed_contig_hits)

        self.assertEqual([1, 1], end_reaches)
        self.assertEqual([2, 1], end_sums)
        self.assertEqual([[1, 1], [0, 1]], end_reached_matrix)
        self.assertEqual([['Contig_1', 450, 750, 'Primer_1', '0'], ['Contig_1', 600, 750, 'Primer_2', '1']], intervals)

    def test_two_seed_sequence_reach_both_ends(self):
        ''' Test that two seed sequences that reach both ends are correctly identified '''
        genome_file_fai = 'TestPrimerReachContigEndCalculation/single_contig_1200N.fasta.fai'
        max_distance = 1000
        seed_contig_hits = {'Contig_1': [['Contig_1', 450, 750, 'Primer_1'], ['Contig_1', 600, 750, 'Primer_2']]}

        end_reaches, end_sums, end_reached_matrix, intervals = \
            search_insertion_sites.seed_reach_contig_end_calc(genome_file_fai, max_distance, seed_contig_hits)

        self.assertEqual([2, 0], end_reaches)
        self.assertEqual([2, 2], end_sums)
        self.assertEqual([[1, 1], [1, 1]], end_reached_matrix)
        self.assertEqual([['Contig_1', 450, 750, 'Primer_1', '0'], ['Contig_1', 600, 750, 'Primer_2', '1']], intervals)


class TestFlankingRegion(unittest.TestCase):

    def test_no_max_distance_limit(self):
        ''' Test that the correct evidence level is returned when no max limit is given. '''
        with open('TestFlankingRegion/single_contig/single_contig_multi_hit.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)

        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 0, 'test')

        self.assertEqual(2, flanking_return)

    def test_multiple_hit_single_contig_w_no_overlaps(self):
        ''' Test the handling of multiple hits from both seed sequneces in a pair but no connection between them '''
        with open('TestFlankingRegion/single_contig/single_contig_multi_hit.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/single_contig/single_contig_1200N.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 1, genome_fai_file)

        self.assertEqual(1, flanking_return)

    def test_multiple_hit_single_contig_w_single_overlap(self):
        ''' Test that seed seqeunces from pair with multiple hit on single contig can be connected correctly with correct max distance'''
        with open('TestFlankingRegion/single_contig/single_contig_multi_hit.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/single_contig/single_contig_1200N.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 51, genome_fai_file)

        self.assertEqual('5B', flanking_return)

    def test_multiple_hit_single_contig_w_multiple_overlaps(self):
        ''' Test that multiple seed sequnces from pair on same contig, that can all be connected give the right evidence level '''
        with open('TestFlankingRegion/single_contig/single_contig_multi_hit.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/single_contig/single_contig_1200N.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 101, genome_fai_file)

        self.assertEqual(2, flanking_return)

    def test_multiple_hit_single_contig_w_same_seed_single_overlap_and_mix_pair(self):
        ''' Test that when two seed seqeunces overlap they can still be recognised as connected. '''
        with open('TestFlankingRegion/single_contig/single_contig_multi_hit_same_overlap_n_mix_pair.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/single_contig/single_contig_1200N.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 1, genome_fai_file)

        self.assertEqual('5B', flanking_return)

    def test_multiple_hit_single_contig_w_same_seed_multiple_overlap_and_mix_pair(self):
        ''' Test that multiple seed seqeunces on the same contig can be connected even when two seeds from a pair overlap '''
        with open('TestFlankingRegion/single_contig/single_contig_multi_hit_same_overlap_n_mix_pair.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/single_contig/single_contig_1200N.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 21, genome_fai_file)

        self.assertEqual(2, flanking_return)

    def test_single_hits_multiple_contigs_overlap_across_contig(self):
        ''' Test that seed sequences from a pair can be connected across the gap between contigs, if given appropriate max distance '''
        with open('TestFlankingRegion/double_contig/multi_contig_single_pair_hit_across_contig.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/double_contig/double_contig.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 101, genome_fai_file)

        self.assertEqual('4B', flanking_return)

    def test_single_hits_multiple_contigs_no_end_reaced(self):
        ''' Test that given a too little max distance with two seed sequences on separate contigs the correct evidence level is returned  '''
        with open('TestFlankingRegion/double_contig/multi_contig_single_pair_hit_across_contig.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/double_contig/double_contig.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 100, genome_fai_file)

        self.assertEqual(1, flanking_return)

    def test_single_hits_multiple_contigs_all_ends_reaced(self):
        ''' Test that given a too large max distance with two seed sequences on separate contigs returns the correct evidence level '''
        with open('TestFlankingRegion/double_contig/multi_contig_single_pair_hit_across_contig.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/double_contig/double_contig.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 10000, genome_fai_file)

        self.assertEqual(2, flanking_return)

    def test_single_hits_multiple_contigs_two_and_one_ends_reaced(self):
        ''' Test the outcome with two seed sequences on seperate contigs, with max distance reaching one end or two ends of contig, depending on seed sequence'''
        with open('TestFlankingRegion/double_contig/multi_contig_single_pair_hit_across_contig.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/double_contig/double_contig.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 300, genome_fai_file)

        self.assertEqual(2, flanking_return)

    def test_single_hits_single_contig_w_multiple_overlaps_w_ends(self):
        ''' Test that two SS that overlap with one reaching an end gives the correct evidence level '''
        with open('TestFlankingRegion/single_contig/Single_contig_multi_hit_two_ss.json', 'r') as seed_hit_json:
            seed_hit_dict = json.load(seed_hit_json)
        genome_fai_file = 'TestFlankingRegion/single_contig/single_contig_1200N.fasta.fai'
        flanking_return = search_insertion_sites.examine_flanking_regions(seed_hit_dict, 351, genome_fai_file)

        self.assertEqual(2, flanking_return)


class TestWriteBedFromseeds(unittest.TestCase):

    def test_writing_two_seeds(self):
        ''' Test function for writing a bed file from a set of seeds '''
        list_of_seeds = [['Contig_1', '500', '600', 'Primer_1'], ['Contig_1', '600', '700', 'Primer_2']]
        bed_file_name = 'TestWriteBedFromPrimers/unit_test_file.bed'

        search_insertion_sites.write_bed_from_list_of_seeds(list_of_seeds, bed_file_name)

        expected_bed_file = 'TestWriteBedFromPrimers/expected_bed_file.bed'

        with open(expected_bed_file, 'r') as expected:
            with open(bed_file_name, 'r') as result:
                self.assertEqual(expected.readlines(), result.readlines())

        os.remove(bed_file_name)


class TestBedMergeHandling(unittest.TestCase):
    def test_single_connection_of_seed_sequences_exclude_seeds(self):
        ''' Test the merge of two seed sequences on single contig and exclude seeds '''
        blast_hit_beds = ['TestBedMergeHandling/Contig_1~~simple_connect.bed']
        include_seeds = False
        exclude_seed_list = ['TestBedMergeHandling/Contig_1~~simple_connect_exclude_primers_file.bed']
        max_seed_dist = 101
        seed_evidence = {'simple_connect': '5A'}

        merged_bed_files, seed_evidence = search_insertion_sites.bed_merge_handling(blast_hit_beds,
                                                  include_seeds,
                                                  exclude_seed_list,
                                                  max_seed_dist,
                                                  seed_evidence)

        self.assertEqual('5B', seed_evidence['simple_connect'])

        with open(merged_bed_files[0], 'r') as result:
            self.assertEqual(['Contig_1\t600\t700\tsimple_connect_1,simple_connect_2\t2\n'], result.readlines())

        os.remove('TestBedMergeHandling/Contig_1~~simple_connect_merged.bed')

    def test_single_connection_of_seed_sequences_include_seeds(self):
        ''' Test merge of two seed sequences on single contig with inclution of seed seqeunces '''
        blast_hit_beds = ['TestBedMergeHandling/Contig_1~~simple_connect.bed']
        include_seeds = True
        exclude_seed_list = ['TestBedMergeHandling/Contig_1~~simple_connect_exclude_primers_file.bed']
        max_seed_dist = 101
        seed_evidence = {'simple_connect': '5B'}

        merged_bed_files, seed_evidence = search_insertion_sites.bed_merge_handling(blast_hit_beds,
                                                  include_seeds,
                                                  exclude_seed_list,
                                                  max_seed_dist,
                                                  seed_evidence)

        self.assertEqual('5B', seed_evidence['simple_connect'])

        with open(merged_bed_files[0], 'r') as result:
            self.assertEqual(['Contig_1\t500\t800\tsimple_connect_1,simple_connect_2\t2\n'], result.readlines())

        os.remove('TestBedMergeHandling/Contig_1~~simple_connect_merged.bed')

    def test_overlap_connection_of_seed_sequences_exclude_seeds(self):
        ''' Test merge of seed sequences that are already overlapping and exclude the seeds '''
        blast_hit_beds = ['TestBedMergeHandling/Contig_1~~overlap_connect.bed']
        include_seeds = False
        exclude_seed_list = ['TestBedMergeHandling/Contig_1~~overlap_connect_exclude_primers_file.bed']
        max_seed_dist = 101
        seed_evidence = {'overlap_connect': '5B'}

        merged_bed_files, seed_evidence = search_insertion_sites.bed_merge_handling(blast_hit_beds,
                                                  include_seeds,
                                                  exclude_seed_list,
                                                  max_seed_dist,
                                                  seed_evidence)

        self.assertEqual(3, seed_evidence['overlap_connect'])

    def test_overlap_connection_of_seed_sequences_include_seeds(self):
        ''' test the merge of already overlapping seed sequence and include the seeds '''
        blast_hit_beds = ['TestBedMergeHandling/Contig_1~~overlap_connect.bed']
        include_seeds = True
        exclude_seed_list = ['TestBedMergeHandling/Contig_1~~overlap_connect_exclude_primers_file.bed']
        max_seed_dist = 101
        seed_evidence = {'overlap_connect': 7}

        merged_bed_files, seed_evidence = search_insertion_sites.bed_merge_handling(blast_hit_beds,
                                                  include_seeds,
                                                  exclude_seed_list,
                                                  max_seed_dist,
                                                  seed_evidence)

        self.assertEqual('5B', seed_evidence['overlap_connect'])

        with open(merged_bed_files[0], 'r') as result:
            self.assertEqual(['Contig_1\t500\t800\toverlap_connect_1,overlap_connect_2\t2\n'], result.readlines())

        os.remove('TestBedMergeHandling/Contig_1~~overlap_connect_merged.bed')

    def test_merge_with_seed_on_contig_edge_exclude_seed(self):
        ''' Test the merge of seeds where one seed falls on the edge of a contig '''
        blast_hit_beds = ['TestBedMergeHandling/double_contig~~primer_edge_placement.bed']
        include_seeds = False
        exclude_seed_list = ['TestBedMergeHandling/double_contig~~primer_edge_primers.bed']
        max_seed_dist = 101
        seed_evidence = {'primer_edge_placement': 3}

        merged_bed_files, seed_evidence = search_insertion_sites.bed_merge_handling(blast_hit_beds,
                                                                                      include_seeds,
                                                                                      exclude_seed_list,
                                                                                      max_seed_dist,
                                                                                      seed_evidence)

        self.assertEqual(3, seed_evidence['primer_edge_placement'])

        with open(merged_bed_files[0], 'r') as result:
            self.assertEqual(['Contig_2\t0\t100\tprimer_edge_placement_2\t1\n'], result.readlines())

        os.remove('TestBedMergeHandling/double_contig~~primer_edge_placement_merged.bed')


class TestExtractSeqsNAnnots(unittest.TestCase):
    def test_same_contig_fasta_extraction(self):
        ''' Test extraction of fasta sequence given a fasta file with region of interest on single same contig. '''
        merged_bed_files = ['TestExtractSeqsNAnnots/No_extraction/Single_contig/Single_contig~~Single_contig_primer.bed']
        file_type = 'fasta'
        genome_file = 'TestExtractSeqsNAnnots/No_extraction/Single_contig/Single_contig.fna'
        annotation_file = ''
        tmp_folder = 'TestExtractSeqsNAnnots/No_extraction/Single_contig'
        out_path = 'TestExtractSeqsNAnnots/No_extraction/Single_contig'
        seed_pairs = {'Single_contig_primer': ['Single_contig_primer_1', 'Single_contig_primer_2']}
        seed_evidence = {'Single_contig_primer': '5B'}
        print_seq_out = 'All'

        copyfile(genome_file, genome_file+'_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file+'_original', genome_file)

        self.assertEqual(1, len(annots_pr_interval))
        self.assertEqual(None, annots_pr_interval['Single_contig_primer'])

        self.assertEqual(0, len(break_seed_sequence_seeds))

        self.assertEqual('5B', seed_sequence_evidence['Single_contig_primer'])

        self.assertEqual(10, inter_seed_sequence_dist['Single_contig_primer'])

        with open(os.path.join(out_path, 'Single_contig--Single_contig_primer.fasta')) as extracted_fasta:
            self.assertEqual('TTTTTTTTTT\n', extracted_fasta.readlines()[1])

        os.remove(os.path.join(out_path, 'Single_contig--Single_contig_primer.fasta'))

    def test_extraction_from_gff_no_annotation_extracted(self):
        ''' Test extraction of fasta seqeunce given a gff file, with no annotations in region of interest and on the same contig '''
        merged_bed_files = ['TestExtractSeqsNAnnots/No_extraction/Single_contig/Single_contig~~Single_contig_primer.bed']
        file_type = 'gff'
        genome_file = 'TestExtractSeqsNAnnots/No_extraction/Single_contig/Single_contig.fna'
        annotation_file = 'TestExtractSeqsNAnnots/No_extraction/Single_contig/Single_contig_annotations.gff'
        tmp_folder = 'TestExtractSeqsNAnnots/No_extraction/Single_contig'
        out_path = 'TestExtractSeqsNAnnots/No_extraction/Single_contig'
        seed_pairs = {'Single_contig_primer': ['Single_contig_primer_1', 'Single_contig_primer_2']}
        seed_evidence = {'Single_contig_primer': '5B'}
        print_seq_out = 'All'

        copyfile(genome_file, genome_file + '_original')
        copyfile(annotation_file, annotation_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)
        os.rename(annotation_file + '_original', annotation_file)

        self.assertEqual(1, len(annots_pr_interval))
        self.assertEqual(0, annots_pr_interval['Single_contig_primer'])

        self.assertEqual(0, len(break_seed_sequence_seeds))

        self.assertEqual('5B', seed_sequence_evidence['Single_contig_primer'])

        self.assertEqual(10, inter_seed_sequence_dist['Single_contig_primer'])

        with open(os.path.join(out_path, 'Single_contig--Single_contig_primer.fasta')) as extracted_fasta:
            self.assertEqual('TTTTTTTTTT\n', extracted_fasta.readlines()[1])

        os.remove(os.path.join(out_path, 'Single_contig--Single_contig_primer.fasta'))

    def test_extraction_from_gff_with_annotation_extracted(self):
        ''' Test extraction of fasta and annotations on a single contig given a gff file '''
        merged_bed_files = ['TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contig_extract_annots~~Single_contig_primer.bed']
        file_type = 'gff'
        genome_file = 'TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contig_extract_annots.fna'
        annotation_file = 'TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contigs_extract_annots.gff'
        tmp_folder = 'TestExtractSeqsNAnnots/With_extraction/Single_contig'
        out_path = 'TestExtractSeqsNAnnots/With_extraction/Single_contig'
        seed_pairs = {'Single_contig_primer': ['Single_contig_primer_1', 'Single_contig_primer_2']}
        seed_evidence = {'Single_contig_primer': '5B'}
        print_seq_out = 'All'

        copyfile(genome_file, genome_file + '_original')
        copyfile(annotation_file, annotation_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)
        os.rename(annotation_file + '_original', annotation_file)

        self.assertEqual(1, len(annots_pr_interval))
        self.assertEqual(2, annots_pr_interval['Single_contig_primer'])

        self.assertEqual(0, len(break_seed_sequence_seeds))

        self.assertEqual('5C', seed_sequence_evidence['Single_contig_primer'])

        self.assertEqual(20, inter_seed_sequence_dist['Single_contig_primer'])

        with open(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.fasta')) as extracted_fasta:
            self.assertEqual('NNNNNNTTTTTTTTTTNNNN\n', extracted_fasta.readlines()[1])

        # Test placement of genes gff when coordinates are adjusted.
        with open(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.gff')) as extracted_annotations:
            gff_lines = extracted_annotations.readlines()
            fist_line = gff_lines[2].split('\t')
            second_line = gff_lines[3].split('\t')
            self.assertListEqual(['1', '3'], fist_line[3:5])
            self.assertListEqual(['17', '19'], second_line[3:5])

            self.assertEqual('NNNNNNTTTTTTTTTTNNNN\n', gff_lines[6])

        os.remove(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.fasta'))
        os.remove(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.gff'))

    def test_cross_contig_fasta_extraction(self):
        ''' Test extraction of fasta sequence across contigs given a fasta file '''
        merged_bed_files = ['TestExtractSeqsNAnnots/No_extraction/Cross_contig/Multi_contig~~Multi_contig_primer.bed']
        file_type = 'fasta'
        genome_file = 'TestExtractSeqsNAnnots/No_extraction/Cross_contig/Multi_contig.fna'
        annotation_file = ''
        tmp_folder = 'TestExtractSeqsNAnnots/No_extraction/Cross_contig'
        out_path = 'TestExtractSeqsNAnnots/No_extraction/Cross_contig'
        seed_pairs = {'Multi_contig_primer': ['Multi_contig_primer_1', 'Multi_contig_primer_2']}
        seed_evidence = {'Multi_contig_primer': '4B'}
        print_seq_out = 'All'

        copyfile(genome_file, genome_file+'_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file+'_original', genome_file)

        self.assertEqual(1, len(annots_pr_interval))
        self.assertEqual(None, annots_pr_interval['Multi_contig_primer'])

        self.assertEqual(2, len(break_seed_sequence_seeds))
        expected_seeds_returned = {'Multi_contig_primer_1_break': ['Multi_contig_primer_1', 'break'], 'Multi_contig_primer_2_break': ['Multi_contig_primer_2', 'break']}
        self.assertEqual(expected_seeds_returned, break_seed_sequence_seeds)

        self.assertEqual('4B', seed_sequence_evidence['Multi_contig_primer'])

        self.assertEqual(10, inter_seed_sequence_dist['Multi_contig_primer_1_break'])
        self.assertEqual(10, inter_seed_sequence_dist['Multi_contig_primer_2_break'])

        with open(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_1_break.fasta')) as extracted_fasta:
            self.assertEqual('AAAAAAAAAA\n', extracted_fasta.readlines()[1])

        with open(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_2_break.fasta')) as extracted_fasta:
            self.assertEqual('AAAAAAAAAA\n', extracted_fasta.readlines()[1])

        os.remove(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_1_break.fasta'))
        os.remove(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_2_break.fasta'))

    def test_extraction_from_gff_across_contigs_no_annotation_extracted(self):
        ''' Test extraction of fasta sequence across contigs given a gff file with no annotations in the region to be extracted '''
        merged_bed_files = ['TestExtractSeqsNAnnots/No_extraction/Cross_contig/Multi_contig~~Multi_contig_primer.bed']
        file_type = 'gff'
        genome_file = 'TestExtractSeqsNAnnots/No_extraction/Cross_contig/Multi_contig.fna'
        annotation_file = 'TestExtractSeqsNAnnots/No_extraction/Cross_contig/Multi_contig_annotations.gff'
        tmp_folder = 'TestExtractSeqsNAnnots/No_extraction/Cross_contig'
        out_path = 'TestExtractSeqsNAnnots/No_extraction/Cross_contig'
        seed_pairs = {'Multi_contig_primer': ['Multi_contig_primer_1', 'Multi_contig_primer_2']}
        seed_evidence = {'Multi_contig_primer': '4B'}
        print_seq_out = 'All'

        copyfile(genome_file, genome_file + '_original')
        copyfile(annotation_file, annotation_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)
        os.rename(annotation_file + '_original', annotation_file)

        self.assertEqual(3, len(annots_pr_interval))
        expected_annots = {'Multi_contig_primer': None, 'Multi_contig_primer_1_break': 0, 'Multi_contig_primer_2_break': 0}
        self.assertEqual(expected_annots, annots_pr_interval)

        self.assertEqual(2, len(break_seed_sequence_seeds))

        self.assertEqual('4B', seed_sequence_evidence['Multi_contig_primer'])

        self.assertEqual(10, inter_seed_sequence_dist['Multi_contig_primer_1_break'])
        self.assertEqual(10, inter_seed_sequence_dist['Multi_contig_primer_2_break'])

        with open(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_1_break.fasta')) as extracted_fasta:
            self.assertEqual('AAAAAAAAAA\n', extracted_fasta.readlines()[1])

        with open(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_2_break.fasta')) as extracted_fasta:
            self.assertEqual('AAAAAAAAAA\n', extracted_fasta.readlines()[1])

        os.remove(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_1_break.fasta'))
        os.remove(os.path.join(out_path, 'Multi_contig--Multi_contig_primer_2_break.fasta'))

    def test_extraction_from_gff_across_contigs_with_annotation_extracted(self):
        ''' Test extraction of annotations and fasta sequences given a gff file, with annotations in the region to be extracted across contigs. '''
        merged_bed_files = ['TestExtractSeqsNAnnots/With_extraction/Cross_contig/Multi_contig_extraction~~Multi_contig_extraction_primer.bed']
        file_type = 'gff'
        genome_file = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig/Multi_contig_extraction.fna'
        annotation_file = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig/Multi_contig_extraction_annotations.gff'
        tmp_folder = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig'
        out_path = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig'
        seed_pairs = {'Multi_contig_extraction_primer': ['Multi_contig_extraction_primer_1', 'Multi_contig_extraction_primer_2']}
        seed_evidence = {'Multi_contig_extraction_primer': '4B'}
        print_seq_out = 'All'

        copyfile(genome_file, genome_file + '_original')
        copyfile(annotation_file, annotation_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)
        os.rename(annotation_file + '_original', annotation_file)

        self.assertEqual(3, len(annots_pr_interval))
        expected_annots = {'Multi_contig_extraction_primer': None, 'Multi_contig_extraction_primer_1_break': 1, 'Multi_contig_extraction_primer_2_break': 1}
        self.assertEqual(expected_annots, annots_pr_interval)

        self.assertEqual(2, len(break_seed_sequence_seeds))

        self.assertEqual('4C', seed_sequence_evidence['Multi_contig_extraction_primer'])

        self.assertEqual(14, inter_seed_sequence_dist['Multi_contig_extraction_primer_1_break'])
        self.assertEqual(14, inter_seed_sequence_dist['Multi_contig_extraction_primer_2_break'])

        with open(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_1_break.fasta')) as extracted_fasta:
            self.assertEqual('AAAAAAAAAANNNN\n', extracted_fasta.readlines()[1])

        with open(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_2_break.fasta')) as extracted_fasta:
            self.assertEqual('NNNNAAAAAAAAAA\n', extracted_fasta.readlines()[1])

        os.remove(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_1_break.fasta'))
        os.remove(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_2_break.fasta'))
        os.remove(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_1_break.gff'))
        os.remove(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_2_break.gff'))

    def test_extraction_when_nothing_is_outputted_fasta(self):
        ''' Test extraction of fasta sequence given a fasta file with region of interest on single same contig.
         However, nothing should be printed to output files'''
        merged_bed_files = [
            'TestExtractSeqsNAnnots/No_extraction/Single_contig/Single_contig~~Single_contig_primer.bed']
        file_type = 'fasta'
        genome_file = 'TestExtractSeqsNAnnots/No_extraction/Single_contig/Single_contig.fna'
        annotation_file = ''
        tmp_folder = 'TestExtractSeqsNAnnots/No_extraction/Single_contig'
        out_path = 'TestExtractSeqsNAnnots/No_extraction/Single_contig'
        seed_pairs = {'Single_contig_primer': ['Single_contig_primer_1', 'Single_contig_primer_2']}
        seed_evidence = {'Single_contig_primer': '5B'}
        print_seq_out = 'None'

        copyfile(genome_file, genome_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)

        self.assertEqual(1, len(annots_pr_interval))
        self.assertEqual(None, annots_pr_interval['Single_contig_primer'])

        self.assertEqual(0, len(break_seed_sequence_seeds))

        self.assertEqual('5B', seed_sequence_evidence['Single_contig_primer'])

        self.assertEqual(10, inter_seed_sequence_dist['Single_contig_primer'])

        with self.assertRaises(FileNotFoundError):
            with open(os.path.join(out_path, 'Single_contig--Single_contig_primer.fasta')) as extracted_fasta:
                pass

    def test_extraction_when_nothing_is_outputted_gff(self):
        merged_bed_files = [
            'TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contig_extract_annots~~Single_contig_primer.bed']
        file_type = 'gff'
        genome_file = 'TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contig_extract_annots.fna'
        annotation_file = 'TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contigs_extract_annots.gff'
        tmp_folder = 'TestExtractSeqsNAnnots/With_extraction/Single_contig'
        out_path = 'TestExtractSeqsNAnnots/With_extraction/Single_contig'
        seed_pairs = {'Single_contig_primer': ['Single_contig_primer_1', 'Single_contig_primer_2']}
        seed_evidence = {'Single_contig_primer': '5B'}
        print_seq_out = 'None'

        copyfile(genome_file, genome_file + '_original')
        copyfile(annotation_file, annotation_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)
        os.rename(annotation_file + '_original', annotation_file)

        self.assertEqual(1, len(annots_pr_interval))
        self.assertEqual(2, annots_pr_interval['Single_contig_primer'])

        self.assertEqual(0, len(break_seed_sequence_seeds))

        self.assertEqual('5C', seed_sequence_evidence['Single_contig_primer'])

        self.assertEqual(20, inter_seed_sequence_dist['Single_contig_primer'])

        with self.assertRaises(FileNotFoundError):
            with open(os.path.join(out_path,
                                   'Single_contig_extract_annots--Single_contig_primer.fasta')) as extracted_fasta:
                pass

        with self.assertRaises(FileNotFoundError):
            with open(os.path.join(out_path,
                                   'Single_contig_extract_annots--Single_contig_primer.gff')) as extracted_annotations:
                pass

    def test_no_break_extraction_from_gff_across_contigs_with_annotation(self):
        ''' Test extraction of annotations and fasta sequences given a gff file, with annotations in the region to be extracted across contigs. '''
        merged_bed_files = ['TestExtractSeqsNAnnots/With_extraction/Cross_contig/Multi_contig_extraction~~Multi_contig_extraction_primer.bed']
        file_type = 'gff'
        genome_file = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig/Multi_contig_extraction.fna'
        annotation_file = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig/Multi_contig_extraction_annotations.gff'
        tmp_folder = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig'
        out_path = 'TestExtractSeqsNAnnots/With_extraction/Cross_contig'
        seed_pairs = {'Multi_contig_extraction_primer': ['Multi_contig_extraction_primer_1', 'Multi_contig_extraction_primer_2']}
        seed_evidence = {'Multi_contig_extraction_primer': '4B'}
        print_seq_out = 'output'

        copyfile(genome_file, genome_file + '_original')
        copyfile(annotation_file, annotation_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)
        os.rename(annotation_file + '_original', annotation_file)

        self.assertEqual(3, len(annots_pr_interval))
        expected_annots = {'Multi_contig_extraction_primer': None, 'Multi_contig_extraction_primer_1_break': 1, 'Multi_contig_extraction_primer_2_break': 1}
        self.assertEqual(expected_annots, annots_pr_interval)

        self.assertEqual(2, len(break_seed_sequence_seeds))

        self.assertEqual('4C', seed_sequence_evidence['Multi_contig_extraction_primer'])

        self.assertEqual(14, inter_seed_sequence_dist['Multi_contig_extraction_primer_1_break'])
        self.assertEqual(14, inter_seed_sequence_dist['Multi_contig_extraction_primer_2_break'])

        with self.assertRaises(FileNotFoundError):
            with open(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_1_break.fasta')) as extracted_fasta:
                pass

        with self.assertRaises(FileNotFoundError):
            with open(os.path.join(out_path, 'Multi_contig_extraction--Multi_contig_extraction_primer_2_break.fasta')) as extracted_fasta:
                pass

        with self.assertRaises(FileNotFoundError):
            with open(os.path.join(out_path,
                                   'Multi_contig_extraction - -Multi_contig_extraction_primer_1_break.fasta')) as extracted_annotations:
                pass

        with self.assertRaises(FileNotFoundError):
            with open(os.path.join(out_path,
                                   'Multi_contig_extraction - -Multi_contig_extraction_primer_2_break.fasta')) as extracted_annotations:
                pass

    def test_extraction_from_gff_with_annotation_extracted_no_break_outputs(self):
        ''' Test extraction of fasta and annotations on a single contig given a gff file '''
        merged_bed_files = ['TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contig_extract_annots~~Single_contig_primer.bed']
        file_type = 'gff'
        genome_file = 'TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contig_extract_annots.fna'
        annotation_file = 'TestExtractSeqsNAnnots/With_extraction/Single_contig/Single_contigs_extract_annots.gff'
        tmp_folder = 'TestExtractSeqsNAnnots/With_extraction/Single_contig'
        out_path = 'TestExtractSeqsNAnnots/With_extraction/Single_contig'
        seed_pairs = {'Single_contig_primer': ['Single_contig_primer_1', 'Single_contig_primer_2']}
        seed_evidence = {'Single_contig_primer': '5B'}
        print_seq_out = 'output'

        copyfile(genome_file, genome_file + '_original')
        copyfile(annotation_file, annotation_file + '_original')

        annots_pr_interval, break_seed_sequence_seeds, seed_sequence_evidence, inter_seed_sequence_dist = \
            search_insertion_sites.extract_seqs_n_annots(merged_bed_files, file_type, genome_file,
                                                         annotation_file, tmp_folder, out_path,
                                                         seed_pairs, seed_evidence, print_seq_out)

        os.rename(genome_file + '_original', genome_file)
        os.rename(annotation_file + '_original', annotation_file)

        self.assertEqual(1, len(annots_pr_interval))
        self.assertEqual(2, annots_pr_interval['Single_contig_primer'])

        self.assertEqual(0, len(break_seed_sequence_seeds))

        self.assertEqual('5C', seed_sequence_evidence['Single_contig_primer'])

        self.assertEqual(20, inter_seed_sequence_dist['Single_contig_primer'])

        with open(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.fasta')) as extracted_fasta:
            self.assertEqual('NNNNNNTTTTTTTTTTNNNN\n', extracted_fasta.readlines()[1])

        # Test placement of genes gff when coordinates are adjusted.
        with open(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.gff')) as extracted_annotations:
            gff_lines = extracted_annotations.readlines()
            fist_line = gff_lines[2].split('\t')
            second_line = gff_lines[3].split('\t')
            self.assertListEqual(['1', '3'], fist_line[3:5])
            self.assertListEqual(['17', '19'], second_line[3:5])

            self.assertEqual('NNNNNNTTTTTTTTTTNNNN\n', gff_lines[6])

        os.remove(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.fasta'))
        os.remove(os.path.join(out_path, 'Single_contig_extract_annots--Single_contig_primer.gff'))


class TestPartitionOutputs(unittest.TestCase):
    # def setUp(self):
    @classmethod
    def setUpClass(cls):
        ''' Construct the names for different mock output files in fasta and gff format '''
        genomes = ['W4rpi', 'lxO0f', 'UM1Dz', '4gDNy', '3PSQZ', 'JltLP']
        cls.seeds = ['ycFQk', '8VNvY', 'Tl04Z', '4EBZ0', 'qngws', 'J08Tv', 'A', 'AA', 'AAA']

        cls.random_fastas = dict.fromkeys(cls.seeds)
        cls.random_gffs = dict.fromkeys(cls.seeds)
        cls.random_breaks = dict.fromkeys(cls.seeds)

        for seed in cls.seeds:
            cls.random_fastas[seed] = []
            cls.random_gffs[seed] = []
            cls.random_breaks[seed] = []
            for genome in genomes:
                cls.random_fastas[seed].append(f'{genome}--{seed}.fasta')
                cls.random_gffs[seed].append(f'{genome}--{seed}.gff')
                cls.random_breaks[seed].append(f'{genome}--{seed}_1_break.fasta')
                cls.random_breaks[seed].append(f'{genome}--{seed}_2_break.fasta')

        os.remove('TestPartitionOutputs/file.txt')

        # Construct mock logger
        cls.logger = logging.getLogger('test_logger.log')
        cls.logger.setLevel(logging.INFO)

    def tearDown(self):
        ''' Class to remove the mock output files and folders '''
        for seed in self.seeds:
            work_dir = os.path.join('TestPartitionOutputs', seed)
            for file in os.listdir(work_dir):
                os.remove(os.path.join(work_dir, file))

            os.rmdir(work_dir)

    def test_partitioning_of_output_fasta_files(self):
        ''' Test that the function for partitioning fasta outputs into separate folder is working as intended '''
        unittest_data_dir = 'TestPartitionOutputs'
        # Construct files to be sorted
        for seed in self.random_fastas:
            for file_name in self.random_fastas[seed]:
                with open(os.path.join(unittest_data_dir, file_name), 'w'):
                    pass

        # remove .DS_Store if on Mac....
        try:
            os.remove(os.path.join(unittest_data_dir, '.DS_Store'))
        except FileNotFoundError:
            pass

        # Run test to wrangle outputs
        wrangle_outputs.partition_outputs(self.seeds, unittest_data_dir, self.logger)

        # Check that all files are in their expected place
        file_presence = []
        for seed in self.seeds:
            for file in self.random_fastas[seed]:

                file = os.path.join(os.path.join(unittest_data_dir, seed), file.replace('--', '-'))
                file_presence.append(os.path.isfile(file))

        # Test itself
        self.assertTrue(all(file_presence))

    def test_partitioning_of_output_gff_files(self):
        ''' Test that the function for partitioning gff outputs into separate folder is working as intended '''
        unittest_data_dir = 'TestPartitionOutputs'
        # Construct files to be sorted
        for seed in self.random_gffs:
            for file_name in self.random_gffs[seed]:
                with open(os.path.join(unittest_data_dir, file_name), 'w'):
                    pass

        # remove .DS_Store if on Mac....
        try:
            os.remove(os.path.join(unittest_data_dir, '.DS_Store'))
        except FileNotFoundError:
            pass

        # Run test to wrangle outputs
        wrangle_outputs.partition_outputs(self.seeds, unittest_data_dir, self.logger)

        # Check that all files are in their expected place
        file_presence = []
        for seed in self.seeds:
            for file in self.random_gffs[seed]:
                file = os.path.join(os.path.join(unittest_data_dir, seed), file.replace('--', '-'))
                file_presence.append(os.path.isfile(file))

        # Test itself
        self.assertTrue(all(file_presence))

    def test_partitioning_of_output_break_files(self):
        ''' Test that the function for partitioning break outputs into separate folder is working as intended '''
        unittest_data_dir = 'TestPartitionOutputs'
        # Construct files to be sorted
        for seed in self.random_breaks:
            for file_name in self.random_breaks[seed]:
                with open(os.path.join(unittest_data_dir, file_name), 'w'):
                    pass

        # remove .DS_Store if on Mac....
        try:
            os.remove(os.path.join(unittest_data_dir, '.DS_Store'))
        except FileNotFoundError:
            pass

        # Run test to wrangle outputs
        wrangle_outputs.partition_outputs(self.seeds, unittest_data_dir, self.logger)

        # Check that all files are in their expected place
        file_presence = []
        for seed in self.seeds:
            for file in self.random_breaks[seed]:
                file = os.path.join(os.path.join(unittest_data_dir, seed), file.replace('--', '-'))
                file_presence.append(os.path.isfile(file))

        # Test itself
        self.assertTrue(all(file_presence))

    @classmethod
    def tearDownClass(cls):
        ''' Write a placeholder file in the folder until next test '''
        with open('TestPartitionOutputs/file.txt', 'w'):
            pass


class TestWritingOutputFiles(unittest.TestCase):

    def test_writing_seed_pairs(self):
        ''' Test the function for writing the file on how seed sequences were paired '''
        seed_pairs = {'primer_1': ['primer_1_1', 'primer_1_2'],
                        'primer_2': ['primer_2_1', 'primer_2_2']}
        wrangle_outputs.write_paired_seeds(seed_pairs, 'TestWritingOutputFiles')

        with open('TestWritingOutputFiles/primer_pairing.expected', 'r') as expected:
            with open('TestWritingOutputFiles/seed_pairing.tsv') as result:
                self.assertEqual(expected.readlines(), result.readlines())

        os.remove('TestWritingOutputFiles/seed_pairing.tsv')

    def test_writing_write_seed_hit_matrix(self):
        ''' Test the function for writing output for the number of hits by each seed sequence in the genome '''
        master_seed_hits = {'genome_1': {'genome': 'genome_1',
                                           'primer_1': 2,
                                           'primer_2': 4},
                              'genome_2': {'genome': 'genome_2',
                                           'primer_1': 2,
                                           'primer_2': 2}}
        seed_pairs = {'primer_1': ['primer_1_1', 'primer_1_2'],
                        'primer_2': ['primer_2_1', 'primer_2_2']}
        out_path = 'TestWritingOutputFiles'

        write_output_csv.write_seed_hit_matrix(master_seed_hits, seed_pairs, out_path)

        with open('TestWritingOutputFiles/contig_hit_matrix.expected', 'r') as expected:
            with open('TestWritingOutputFiles/contig_hit_matrix.csv') as result:
                self.assertEqual(expected.readlines(), result.readlines())

        os.remove('TestWritingOutputFiles/contig_hit_matrix.csv')

    def test_writing_annotation_matrix(self):
        ''' Test the function for writing output for the number of annotations found in the region of interest '''
        master_annotation_hits = {'genome_1': {'genome': 'genome_1',
                                           'primer_1': 2,
                                           'primer_2': 4},
                                  'genome_2': {'genome': 'genome_2',
                                               'primer_1': 2,
                                               'primer_2': 2}}
        seed_pairs = {'primer_1': ['primer_1_1', 'primer_1_2'],
                        'primer_2': ['primer_2_1', 'primer_2_2']}
        out_path = 'TestWritingOutputFiles'

        write_output_csv.write_annotation_num_matrix(master_annotation_hits, seed_pairs, out_path)

        with open('TestWritingOutputFiles/contig_hit_matrix.expected', 'r') as expected:
            with open('TestWritingOutputFiles/annotation_num_matrix.csv') as result:
                self.assertEqual(expected.readlines(), result.readlines())

        os.remove('TestWritingOutputFiles/annotation_num_matrix.csv')

    def test_writing_seed_evidence(self):
        ''' Test the function for writing output for evidence levels for each seed sequence pair '''
        master_seed_evidence = {'genome_1': {'genome': 'genome_1',
                                           'primer_1': 2,
                                           'primer_2': 4},
                                  'genome_2': {'genome': 'genome_2',
                                               'primer_1': 2,
                                               'primer_2': 2}}
        seed_pairs = {'primer_1': ['primer_1_1', 'primer_1_2'],
                        'primer_2': ['primer_2_1', 'primer_2_2']}
        out_path = 'TestWritingOutputFiles'

        write_output_csv.write_seed_hit_evidence(master_seed_evidence, seed_pairs, out_path)

        with open('TestWritingOutputFiles/contig_hit_matrix.expected', 'r') as expected:
            with open('TestWritingOutputFiles/master_seed_evidence.csv') as result:
                self.assertEqual(expected.readlines(), result.readlines())

        os.remove('TestWritingOutputFiles/master_seed_evidence.csv')

    def test_writing_inter_seed_distance(self):
        ''' Test the function for writing output for the distance between connected seed sequences '''
        master_inter_seed_dist = {'genome_1': {'genome': 'genome_1',
                                               'primer_1': 2,
                                               'primer_2': 4},
                                  'genome_2': {'genome': 'genome_2',
                                               'primer_1': 2,
                                               'primer_2': 2}}
        seed_pairs = {'primer_1': ['primer_1_1', 'primer_1_2'],
                        'primer_2': ['primer_2_1', 'primer_2_2']}
        out_path = 'TestWritingOutputFiles'

        write_output_csv.write_inter_seed_dist(master_inter_seed_dist, seed_pairs, out_path)

        with open('TestWritingOutputFiles/contig_hit_matrix.expected', 'r') as expected:
            with open('TestWritingOutputFiles/inter_seed_distance.csv') as result:
                self.assertEqual(expected.readlines(), result.readlines())

        os.remove('TestWritingOutputFiles/inter_seed_distance.csv')


if __name__ == '__main__':
    unittest.main()

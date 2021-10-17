'''
Unit tests for Magphi.

Usage: python -m unittest -v Magphi_test
'''

import unittest
import os

from Magphi import commandline_interface
from Magphi import check_inputs
from Magphi import split_gff_file
from Magphi import primer_handling

from io import StringIO
#pylint: disable=no-name-in-module
from Magphi.__main__ import FastaStats


class TestCommandLineHelpCalls(unittest.TestCase):
    '''Unit test for the commandline interface'''
    def test_no_input(self):
        with self.assertRaises(SystemExit):
            commandline_interface.get_commandline_arguments([], 1)

    def test_single_dash_help(self):
        with self.assertRaises(SystemExit):
            commandline_interface.get_commandline_arguments('-help', 1)

    def test_unrecognised_argument_exit(self):
        with self.assertRaises(SystemExit):
            commandline_interface.get_commandline_arguments(['-p', 'test.file', '-g', 'test.file', '--none'], 1)


class TestFileRecognition(unittest.TestCase):
    # Move to folder with mock input files. First try Github structure, then try pulled repository structure
    try:
        os.chdir('/Magphi/unit_test_data/TestFileRecognition')
    except FileNotFoundError:
        os.chdir('../unit_test_data/TestFileRecognition')

    def test_fasta_recognition(self):
        ''' test the recognition of fasta files '''
        path = 'Fasta_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files)

        self.assertEqual('fasta', file_type)

    def test_none_fasta_recognition(self):
        ''' test that gff files are not recognised as fasta files '''
        path = 'Gff3_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files)

        self.assertEqual(None, file_type)

    def test_mixed_gff_and_fasta_recognition(self):
        ''' test that a mix of fasta and gff files results in exiting Magphi with an error '''
        path = 'Mixed_gff_and_fasta'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_fasta(files)

    def test_fasta_and_random_text_recognition(self):
        ''' test that a mix of fasta and random text files results in exiting Magphi with an error '''
        path = 'Mixed_fasta_and_text'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_fasta(files)

    def test_complete_gff_recognition(self):
        ''' test that gff files with an attached genome are recognised correctly '''
        path = 'Gff3_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_gff(files)

        self.assertEqual('gff', file_type)

    def test_none_gff_recognition(self):
        ''' test that fasta files are not recognised as gff '''
        path = 'Fasta_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_gff(files)

        self.assertEqual(None, file_type)

    def test_gff_missing_genome_recognition(self):
        ''' test that gff files without a genomes attached exits with an error '''
        path = 'Gff3_without_genome_attached'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]
        with self.assertRaises(SystemExit):
            check_inputs.check_if_gff(files)

    def test_gff_and_random_text_recognition(self):
        ''' test that a mix of GFF3 and random text files results in exiting Magphi with an error '''
        path = 'Mixed_gff_and_text'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_gff(files)

    def test_not_incompatible_recognition(self):
        ''' test that a text file not being a Fasta or GFF3 files results in an error '''
        files = ['Mixed_miscellaneous_files/Random_text.txt']

        with self.assertRaises(SystemExit):
            check_inputs.check_inputs(files)

    def test_empty_file(self):
        ''' Test that an empty file results in an error '''
        files = ['Mixed_miscellaneous_files/empty_file.txt']

        with self.assertRaises(SystemExit):
            check_inputs.check_inputs(files)


class TestSplittingGff(unittest.TestCase):
    def test_gff_split_single_file(self):
        ''' test the function that splits a gff file into annotations and genome. Assess the number of lines in output '''
        # Grab the test file
        try:
            os.chdir('/Magphi/unit_test_data/TestSplittingGff')
        except FileNotFoundError:
            os.chdir('../TestSplittingGff')
            #unit_test_data/TestSplittingGff

        path = os.getcwd()
        file = os.path.join(path, 'minimized.gff')

        # Split the test file
        genome, annotation = split_gff_file.split_single_gff(file, path)

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


class TestPrimerFunctions(unittest.TestCase):
    try:
        os.chdir('/Magphi/unit_test_data/TestPrimerFunctions')
    except FileNotFoundError:
        os.chdir('../TestPrimerFunctions')

    def test_uneven_primer_number(self):
        with self.assertRaises(SystemExit):
            primer_handling.check_number_of_primers('Uneven_number_primers.txt')

    def test_correct_primer_pairing(self):
        primer_names = ['D_1', 'D_2', 'mutsD_1', 'mutsD_2']
        primer_pairs = primer_handling.construct_pair_primers(primer_names)

        expected_names = {'D': ['D_1', 'D_2'], 'mutsD': ['mutsD_1', 'mutsD_2']}

        # Sort dicts to make them same order
        primer_pairs = [primer_pairs[key].sort() for key in primer_pairs]
        expected_names = [expected_names[key].sort() for key in expected_names]

        self.assertEqual(expected_names , primer_pairs)

    def test_non_matching_primer_names(self):
        primer_names = ['D_2', 'B_1', 'Z_1', 'A_1']
        with self.assertRaises(SystemExit):
            primer_handling.construct_pair_primers(primer_names)

    def test_identical_primer_names(self):
        with self.assertRaises(SystemExit):
            primer_handling.extract_primer_info('Same_name_primers.txt')

# Bioinitio tests
# class TestFastaStats(unittest.TestCase):
#     '''Unit tests for FastaStats'''
#     def do_test(self, input_str, minlen, expected):
#         "Wrapper function for testing FastaStats"
#         result = FastaStats().from_file(StringIO(input_str), minlen)
#         self.assertEqual(expected, result)
#
#     def test_zero_byte_input(self):
#         "Test input containing zero bytes"
#         expected = FastaStats(num_seqs=0,
#                               num_bases=0,
#                               min_len=None,
#                               max_len=None,
#                               average=None)
#         self.do_test('', 0, expected)
#
#     def test_single_newline_input(self):
#         "Test input containing a newline (\n) character"
#         expected = FastaStats(num_seqs=0,
#                               num_bases=0,
#                               min_len=None,
#                               max_len=None,
#                               average=None)
#         self.do_test('\n', 0, expected)
#
#     def test_single_greater_than_input(self):
#         "Test input containing a single greater-than (>) character"
#         expected = FastaStats(num_seqs=1,
#                               num_bases=0,
#                               min_len=0,
#                               max_len=0,
#                               average=0)
#         self.do_test('>', 0, expected)
#
#     def test_one_sequence(self):
#         "Test input containing one sequence"
#         expected = FastaStats(num_seqs=1,
#                               num_bases=5,
#                               min_len=5,
#                               max_len=5,
#                               average=5)
#         self.do_test(">header\nATGC\nA", 0, expected)
#
#     def test_two_sequences(self):
#         "Test input containing two sequences"
#         expected = FastaStats(num_seqs=2,
#                               num_bases=9,
#                               min_len=2,
#                               max_len=7,
#                               average=4)
#         self.do_test(">header1\nATGC\nAGG\n>header2\nTT\n", 0, expected)
#
#     def test_no_header(self):
#         "Test input containing sequence without preceding header"
#         expected = FastaStats(num_seqs=0,
#                               num_bases=0,
#                               min_len=None,
#                               max_len=None,
#                               average=None)
#         self.do_test("no header\n", 0, expected)
#
#     def test_minlen_less_than_all(self):
#         "Test input when --minlen is less than 2 out of 2 sequences"
#         expected = FastaStats(num_seqs=2,
#                               num_bases=9,
#                               min_len=2,
#                               max_len=7,
#                               average=4)
#         self.do_test(">header1\nATGC\nAGG\n>header2\nTT\n", 2, expected)
#
#     def test_minlen_greater_than_one(self):
#         "Test input when --minlen is less than 1 out of 2 sequences"
#         expected = FastaStats(num_seqs=1,
#                               num_bases=7,
#                               min_len=7,
#                               max_len=7,
#                               average=7)
#         self.do_test(">header1\nATGC\nAGG\n>header2\nTT\n", 3, expected)
#
#     def test_minlen_greater_than_all(self):
#         "Test input when --minlen is greater than 2 out of 2 sequences"
#         expected = FastaStats(num_seqs=0,
#                               num_bases=0,
#                               min_len=None,
#                               max_len=None,
#                               average=None)
#         self.do_test(">header1\nATGC\nAGG\n>header2\nTT\n", 8, expected)


if __name__ == '__main__':
    unittest.main()

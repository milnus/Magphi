'''
Unit tests for Magphi.

Usage: python -m unittest -v Magphi_test
'''

import unittest
import os

from Magphi import commandline_interface
from Magphi import check_inputs

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
        ''' Function to test the recognition of differently formatted fasta files '''
        path = 'test_fasta_recognition'
        files = os.listdir('test_fasta_recognition')
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files)

        self.assertEqual('fasta', file_type)

    # def test_mixed_recognition(self):
    #     path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_mixed'
    #     files = os.listdir(path)
    #     files = [os.path.join(path, file) for file in files]
    #
    #     with self.assertRaises(SystemExit):
    #         check_inputs.check_if_fasta(files)
    #
    # def test_none_fasta_recognition(self):
    #     path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_GFF3'
    #     files = os.listdir(path)
    #     files = [os.path.join(path, file) for file in files]
    #
    #     file_type = check_inputs.check_if_fasta(files)
    #
    #     self.assertEqual(None, file_type)
    #
    # def test_complete_gff_recognition(self):
    #     path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_GFF3'
    #     files = os.listdir(path)
    #     files = [os.path.join(path, file) for file in files]
    #
    #     file_type = check_inputs.check_if_gff(files)
    #
    #     self.assertEqual('gff', file_type)
    #
    # def test_gff_missing_genome_recognition(self):
    #     path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_GFF3_wo_genome'
    #     files = os.listdir(path)
    #     files = [os.path.join(path, file) for file in files]
    #
    #     with self.assertRaises(SystemExit):
    #         check_inputs.check_if_gff(files)
    #
    # def test_not_incompatible_recognition(self):
    #     path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_not_recognised_files'
    #     files = os.listdir(path)
    #     files = [os.path.join(path, file) for file in files]
    #
    #     with self.assertRaises(SystemExit):
    #         check_inputs.check_inputs(files)
    #
    # def test_empty_file(self):
    #     pass

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

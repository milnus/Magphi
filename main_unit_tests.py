import unittest
import os
import check_inputs
import commandline_interface
import split_gff_file
import primer_handling
import search_insertion_sites
import json


# Test for printing help function upon no inputs
class TestCommandLineHelpCalls(unittest.TestCase):
    def test_dashhelp(self):
        with self.assertRaises(SystemExit):
            commandline_interface.get_commandline_arguments('-help')

    def test_no_input(self):
        with self.assertRaises(SystemExit):
            commandline_interface.get_commandline_arguments([])


class TestFileRecognition(unittest.TestCase):
    def test_fasta_recognition(self):
        path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_fasta'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files)

        self.assertEqual('fasta', file_type)

    def test_mixed_recognition(self):
        path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_mixed'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_fasta(files)

    def test_none_fasta_recognition(self):
        path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_GFF3'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_fasta(files)

        self.assertEqual(None, file_type)

    def test_complete_gff_recognition(self):
        path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_GFF3'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        file_type = check_inputs.check_if_gff(files)

        self.assertEqual('gff', file_type)

    def test_gff_missing_genome_recognition(self):
        path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_recognise_GFF3_wo_genome'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_if_gff(files)

    def test_not_incompatible_recognition(self):
        path = '/Users/mjespersen/Documents/Phupa_test_data/Unittest_not_recognised_files'
        files = os.listdir(path)
        files = [os.path.join(path, file) for file in files]

        with self.assertRaises(SystemExit):
            check_inputs.check_inputs(files)

# Test splitting of gff file
class TestSplittingGff(unittest.TestCase):

    def test_gff_split_single_file(self):
        # Grab the test file
        path = '/Users/mjespersen/Documents/Phupa_test_data/Unitest_split_GFF3'
        files = os.listdir(path)
        files = [file for file in files if file != ".DS_Store"]
        file = os.path.join(path, files[0])

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

    def test_uneven_primer_number(self):
        with self.assertRaises(SystemExit):
            primer_handling.check_number_of_primers('/Users/mjespersen/Documents/Phupa_test_data/Unittest_uneven_primer_number/Untitled_primers.txt')

class TestFlankingRegion(unittest.TestCase):

    def test_multiple_hit_single_contig_no_limit(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 0, 'test')

        self.assertEqual(3, flanking_return)

    def test_multiple_hit_single_contig_w_single_overlap(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 201, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_1200N.fasta.fai')

        self.assertEqual(5, flanking_return)

    def test_multiple_hit_single_contig_w_multiple_overlaps(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 300, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_1200N.fasta.fai')

        self.assertEqual(3, flanking_return)

    def test_multiple_hit_single_contig_w_no_overlaps(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 1, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_1200N.fasta.fai')

        self.assertEqual(2, flanking_return)

    def test_multiple_hit_single_contig_w_same_primer_overlap(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit_sampe_overlap.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 1, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_1200N.fasta.fai')

        self.assertEqual(2, flanking_return)


    def test_multiple_hit_single_contig_w_same_primer_overlap_and_mix_pair(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit_sampe_overlap_n_mix_pair.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 1, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_1200N.fasta.fai')

        self.assertEqual(5, flanking_return)

    def test_multiple_hit_single_contig_w_same_primer_neighbour_and_mate(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit_same_neighbours.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 200, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_1200N.fasta.fai')

        self.assertEqual(5, flanking_return)

    def test_multiple_hit_single_contig_w_same_primer_neighbour_and_mate_multi_hit(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_multi_hit_same_neighbours.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 400, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/single_contig_1200N.fasta.fai')

        self.assertEqual(3, flanking_return)


    def test_multiple_hits_multiple_contigs_overlap_across_contig(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/multi_contig_multi_hit_overlap_across_contig.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 101, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/double_contig.fasta.fai')

        self.assertEqual(6, flanking_return)

    def test_multiple_hits_multiple_contigs_no_end_reaced(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/multi_contig_multi_hit_overlap_across_contig.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 100, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/double_contig.fasta.fai')

        self.assertEqual(2, flanking_return)

    def test_multiple_hits_multiple_contigs_all_ends_reaced(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/multi_contig_multi_hit_overlap_across_contig.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 10000, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/double_contig.fasta.fai')

        self.assertEqual(3, flanking_return)

    def test_multiple_hits_multiple_contigs_two_and_one_ends_reaced(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/multi_contig_multi_hit_overlap_across_contig.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 300, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/double_contig.fasta.fai')

        self.assertEqual(4, flanking_return)

    # TODO - Write a unittest with only one primer that has hit.

    # TODO - Write a unittest with only one primer reaching one end of a contig
    def test_one_primer_reaching_one_end(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/multi_contig_multi_hit_single_primer_reaching_one_end.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 101, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/double_contig.fasta.fai')

        self.assertEqual(6, flanking_return)

    # TODO - Write a unit test with only one primer reaching both ends of a contig
    def test_one_primer_reaching_two_ends(self):
        with open('/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/multi_contig_multi_hit_single_primer_reaching_one_end.json', 'r') as primer_hit_json:
            primer_hit_dict = json.load(primer_hit_json)
        flanking_return = search_insertion_sites.examine_flanking_regions(primer_hit_dict, 400, '/Users/mjespersen/Documents/Phupa_test_data/Unitest_flanking_regions/double_contig.fasta.fai')

        self.assertEqual(6, flanking_return)

    # TODO - test warning return - if you can figure out how to do it ;-)

if __name__ == '__main__':
    unittest.main()
#!/usr/bin/env bash

# 1. Parse command line arguments
# 2. cd to the test directory
# 3. run tests
# 4. Print summary of successes and failures, exit with 0 if
#    all tests pass, else exit with 1

# Uncomment the line below if you want more debugging information
# about this script.
#set -x

# The name of this test script
this_program_name="Magphi-test.sh"
# The program we want to test (either a full path to an executable, or the name of an executable in $PATH)
test_program=""
# Directory containing the test data files and expected outputs
test_data_dir=""
# Number of failed test cases
num_errors=0
# Total number of tests run
num_tests=0

function show_help {
cat << UsageMessage

${this_program_name}: run integration/regression tests for Magphi 

Usage:
    ${this_program_name} [-h] [-v] -p program -d test_data_dir 

Example:
    ${this_program_name} -p bin/Magphi -d data/tests

-h shows this help message

-v verbose output
UsageMessage
}

# echo an error message $1 and exit with status $2
function exit_with_error {
    printf "${this_program_name}: ERROR: $1\n"
    exit $2
}

# if -v is specified on the command line, print a more verbaose message to stdout
function verbose_message {
    if [ "${verbose}" = true ]; then
        echo "${this_program_name} $1"
    fi
}

# Parse the command line arguments and set the global variables program and test_data_dir 
function parse_args {
    local OPTIND opt

    while getopts "hp:d:v" opt; do
        case "${opt}" in
            h)
                show_help
                exit 0
                ;;
            p)  test_program="${OPTARG}"
                ;;
            d)  test_data_dir="${OPTARG}"
                ;;
            v)  verbose=true
                ;;
        esac
    done

    shift $((OPTIND-1))

    [ "$1" = "--" ] && shift

    if [[ -z ${test_program} ]]; then
        exit_with_error "missing command line argument: -p program, use -h for help" 2
    fi

    if [[ -z ${test_data_dir} ]]; then
        exit_with_error "missing command line argument: -d test_data_dir, use -h for help" 2
    fi
}


# Run a command and check that the output is
# exactly equal the contents of a specified file 
# ARG1: command we want to test as a string
# ARG2: a file path containing the expected output
# ARG3: expected exit status
function test_stdout_exit {
    let num_tests+=1
    output=$(eval $1)
    exit_status=$?
    expected_output_file=$2
    expected_exit_status=$3
    verbose_message "Testing stdout and exit status: $1"
    difference=$(diff <(echo "$output") $expected_output_file)
    if [ -n "$difference" ]; then 
        let num_errors+=1
        echo "Test output failed: $1"
        echo "Actual output:"
        echo "$output"
        expected_output=$(cat $2)
        echo "Expected output:"
        echo "$expected_output"
        echo "Difference:"
        echo "$difference"
    elif [ "$exit_status" -ne "$expected_exit_status" ]; then
        let num_errors+=1
        echo "Test exit status failed: $1"
        echo "Actual exit status: $exit_status"
        echo "Expected exit status: $expected_exit_status"
    fi 
}

# Run a command and check that the output file is
# exactly equal the contents of a specified file
# ARG1: A file returned from program after running
# ARG2: a file path containing the expected output
function test_output_file {
    let num_tests+=1
    output=$1
    expected_output_file=$2
    verbose_message "Testing output file: $1"
    verbose_message "Expected file path: $2"
    difference=$(diff $output $expected_output_file)
    if [ -n "$difference" ]; then
        let num_errors+=1
        echo "Test output failed: $1"
        echo "Actual output:"
        cat $output
        expected_output=$(cat $expected_output_file)
        echo "Expected output:"
        echo "$expected_output"
        echo "Difference:"
        echo "$difference"
    fi
}

# Run a command and check that the exit status is 
# equal to an expected value
# exactly equal the contents of a specified file 
# ARG1: command we want to test as a string
# ARG2: expected exit status
# NB: this is mostly for checking erroneous conditions, where the
# exact output message is not crucial, but the exit status is
# important
function test_exit_status {
    let num_tests+=1
    output=$(eval $1)
    exit_status=$?
    expected_exit_status=$2
    verbose_message "Testing exit status: $1"
    if [ "$exit_status" -ne "$expected_exit_status" ]; then
        let num_errors+=1
        echo "Test exit status failed: $1"
        echo "Actual exit status: $exit_status"
        echo "Expected exit status: $expected_exit_status"
    fi 
}

function call_new_test {
  echo ''
  echo $1
}

# 1. Parse command line arguments.
parse_args $@
# 2. Change to test directory
cd $test_data_dir
# 2. Run tests
## Test commandline exit status
# Test output for no arguments
call_new_test "Test output for no arguments"
test_stdout_exit "$test_program" no_input.expected 2
# Test output for -help argument given
call_new_test "Test output for -help argument given"
test_stdout_exit "$test_program -help" no_input.expected 0
# Test exit status for a bad command line invocation
call_new_test "Test exit status for a bad command line invocation"
test_exit_status "$test_program --this_is_not_a_valid_argument > /dev/null 2>&1" 2

## Check exit status for bad input GFF or FASTA files
# Test exit status when input is mixed fasta and gff
call_new_test "Test exit status when input is mixed fasta and gff"
test_exit_status "$test_program -g test_fasta.fna test_GFF.gff -s empty_file > /dev/null 2>&1" 3
# Test exit status when fasta is mixed with random text file
call_new_test "Test exit status when fasta is mixed with random text file"
test_exit_status "$test_program -g test_fasta.fna random_text.txt -s empty_file > /dev/null 2>&1" 3
# Test exit status when gff is mixed with random text file
call_new_test "Test exit status when gff is mixed with random text file"
test_exit_status "$test_program -g test_GFF.gff random_text.txt -s empty_file > /dev/null 2>&1" 3
# Test when just random text file is given as input
call_new_test "Test when just random text file is given as input"
test_exit_status "$test_program -g random_text.txt -s empty_file > /dev/null 2>&1" 3
# Test when empty file is given as input
call_new_test "Test when empty file is given as input"
test_exit_status "$test_program -g empty_file -s empty_file > /dev/null 2>&1" 3

# M  - Test multiple hits where two can be found to connect - same contig - # TODO - written -
# N  - Test multiple hits where two can be found to connect - across contigs. - # TODO - -



# Run test for evidence level when no seed hits are found
call_new_test "Run test for evidence level when no seed hits are found"
Magphi -g evidence_levels_simple_genome.fasta -s no_primers_match_primers.fasta -o test_out_folder > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv no_primers_match_evidence_levels.expected
rm -r test_out_folder

# Run test for evidence level when a single seed hit is found
call_new_test "Test for evidence level when a single seed hit is found"
Magphi -g evidence_levels_simple_genome.fasta -s single_primer_match_primers.fasta -o test_out_folder > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv one_primer_match_evidence_level.expected
rm -r test_out_folder

# Run test for evidence level when two primers hit multiple times on same contig, but cannot connect
call_new_test "Test for evidence level when two primers hit multiple times on same contig, but cannot connect"
Magphi -g evidence_levels_simple_genome.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 1 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv two_primers_multiple_hits_single_contig_no_connections.expected
rm -r test_out_folder

# Run test for evidence level when two primers hit the same contig multiple times and they can connect multiple ways
call_new_test "Test for evidence level when two primers hit the same contig multiple times and they can connect multiple ways"
Magphi -g evidence_levels_simple_genome.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 1000 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv two_primers_multiple_hits_single_contig_multiple_connections.expected
rm -r test_out_folder

# Run test for evidence level when two primers hit multiple contigs multiple times and no connection between them
call_new_test "Test for evidence level when two primers hit multiple contigs multiple times and no connection between them"
Magphi -g evidence_levels_simple_two_contigs.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 1 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv two_primers_two_contigs_multipe_hits_no_connection.expected
rm -r test_out_folder

# run test for evidence level when two seeds hit multiple contigs multiple times and multiple connections
call_new_test "Test for evidence level when two seeds hit multiple contigs multiple times and multiple connections"
Magphi -g evidence_levels_simple_two_contigs.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 1000 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv two_primers_two_contigs_multi_hit_multi_connect.expected
rm -r test_out_folder

# Run test for evidence level when two seeds overlap and are excluded because primers are deleted
call_new_test "Test for evidence level when two seeds overlap and are excluded because primers are deleted"
Magphi -g evidence_levels_overlap_two_contigs.fasta -s overlap_primers.fasta -o test_out_folder -md 1 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_overlapping_seeds.expected
rm -r test_out_folder

# Run test for evidence level when one seed is on the edge of a contig and is deleted due to being excluded
call_new_test "Test for evidence level when one seed is on the edge of a contig and is deleted due to being excluded"
Magphi -g evidence_levels_overlap_two_contigs.fasta -s contig_edge_primers.fasta -o test_out_folder -md 180 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_contig_edge.expected
rm -r test_out_folder

# Run test for evidence level when only two unique seeds hit but can not connect
call_new_test "Test for evidence level when only two unique seeds hit but can not connect"
Magphi -g two_contigs_two_primers_single_hit.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 1 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_simple_hits_cross_contig_no_connect.expected
rm -r test_out_folder

# Run test for evidence level when only two unique seeds hit with connection but no annotation
call_new_test "Test for evidence level when only two unique seeds hit with connection but no annotation"
Magphi -g two_contigs_two_primers_single_hit.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 70 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_simple_hits_cross_contig_connected/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_levels_simple_hits_cross_contig_connected/inter_primer_distance.csv
test_output_file test_out_folder/two_primers_simple/two_contigs_two_primers_single_hit-two_primers_simple_1_break.fasta evidence_levels_simple_hits_cross_contig_connected/two_contigs_two_primers_single_hit-two_primers_simple_1_break.fasta
test_output_file test_out_folder/two_primers_simple/two_contigs_two_primers_single_hit-two_primers_simple_2_break.fasta evidence_levels_simple_hits_cross_contig_connected/two_contigs_two_primers_single_hit-two_primers_simple_2_break.fasta
rm -r test_out_folder

# Run test for evidence level when only two unique seeds hit with connection and annotation
call_new_test "Test for evidence level when only two unique seeds hit with connection and annotation"
Magphi -g two_contigs_two_primers_single_hit.gff -s two_seeds_w_annotation_between.fasta -o test_out_folder -md 375 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_simple_primers_connect_cross_contig_with_annotations/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_levels_simple_primers_connect_cross_contig_with_annotations/inter_primer_distance.csv
test_output_file test_out_folder/annotation_num_matrix.csv evidence_levels_simple_primers_connect_cross_contig_with_annotations/annotation_num_matrix.csv
test_output_file test_out_folder/annotation_primers/two_contigs_two_primers_single_hit-annotation_primers_1_break.fasta evidence_levels_simple_primers_connect_cross_contig_with_annotations/two_contigs_two_primers_single_hit-annotation_primers_1_break.fasta
test_output_file test_out_folder/annotation_primers/two_contigs_two_primers_single_hit-annotation_primers_2_break.fasta evidence_levels_simple_primers_connect_cross_contig_with_annotations/two_contigs_two_primers_single_hit-annotation_primers_2_break.fasta
test_output_file test_out_folder/annotation_primers/two_contigs_two_primers_single_hit-annotation_primers_1_break.gff evidence_levels_simple_primers_connect_cross_contig_with_annotations/two_contigs_two_primers_single_hit-annotation_primers_1_break.gff
test_output_file test_out_folder/annotation_primers/two_contigs_two_primers_single_hit-annotation_primers_2_break.gff evidence_levels_simple_primers_connect_cross_contig_with_annotations/two_contigs_two_primers_single_hit-annotation_primers_2_break.gff
rm -r test_out_folder

# Run test on a gff with a single contig no connection between primers
call_new_test 'Test of a gff with a single contig no connection between primers'
Magphi -g evidence_levels_single_contig.gff -s two_primers_simple_match_primers.fasta -o test_out_folder -md 1 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_single_contig_no_connection/master_primer_evidence.csv
rm -r test_out_folder

# Run test on a gff with a single contig with connection between primers
call_new_test "Test of a gff with a single contig with connection between primers"
Magphi -g evidence_levels_single_contig.gff -s two_primers_simple_match_primers.fasta -o test_out_folder -md 100 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_single_contig_connection/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_levels_single_contig_connection/inter_primer_distance.csv
test_output_file test_out_folder/two_primers_simple/evidence_levels_single_contig-two_primers_simple.fasta evidence_levels_single_contig_connection/evidence_levels_single_contig-two_primers_simple.fasta
test_output_file test_out_folder/annotation_num_matrix.csv evidence_levels_single_contig_connection/annotation_num_matrix.csv
rm -r test_out_folder

# Run test on a gff with a single contig with connection and annotation between primers - exclude seeds
call_new_test "Test of a gff with a single contig with connection and annotation between primers - exclude seeds"
Magphi -g evidence_levels_single_contig.gff -s two_seeds_w_annotation_between.fasta -o test_out_folder -md 750 -is > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_level_single_contigs_connect_annotations_include_seeds/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_level_single_contigs_connect_annotations_include_seeds/inter_primer_distance.csv
test_output_file test_out_folder/annotation_primers/evidence_levels_single_contig-annotation_primers.gff evidence_level_single_contigs_connect_annotations_include_seeds/evidence_levels_single_contig-annotation_primers.gff
test_output_file test_out_folder/annotation_primers/evidence_levels_single_contig-annotation_primers.fasta evidence_level_single_contigs_connect_annotations_include_seeds/evidence_levels_single_contig-annotation_primers.fasta
test_output_file test_out_folder/annotation_num_matrix.csv evidence_level_single_contigs_connect_annotations_include_seeds/annotation_num_matrix.csv
rm -r test_out_folder

# Test output when seeds hit multiple times and one and only one connection can be made between seeds on same contig
call_new_test "Test output when seeds hit multiple times and one and only one connection can be made between seeds on same contig"
Magphi -g evidence_levels_simple_two_contigs.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 130 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_level_multi_hit_multi_contig_single_same_contig_overlap/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_level_multi_hit_multi_contig_single_same_contig_overlap/inter_primer_distance.csv
test_output_file test_out_folder/contig_hit_matrix.csv evidence_level_multi_hit_multi_contig_single_same_contig_overlap/contig_hit_matrix.csv
test_output_file test_out_folder/evidence_levels_simple_two_contigs-two_primers_simple.fasta evidence_level_multi_hit_multi_contig_single_same_contig_overlap/evidence_levels_simple_two_contigs-two_primers_simple.fasta
rm -r test_out_folder

# Test output when seeds hit multiple times and one and only one connection can be made between seeds on different contig
call_new_test "Test output when seeds hit multiple times and one and only one connection can be made between seeds on different contig"
Magphi -g evidence_levels_simple_two_contigs_cross_conenct.fasta -s two_primers_simple_match_primers.fasta -o test_out_folder -md 180 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_level_multi_seed_hit_cross_contig_connect/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_level_multi_seed_hit_cross_contig_connect/inter_primer_distance.csv
test_output_file test_out_folder/contig_hit_matrix.csv evidence_level_multi_seed_hit_cross_contig_connect/contig_hit_matrix.csv
test_output_file test_out_folder/evidence_levels_simple_two_contigs-two_primers_simple.fasta evidence_level_multi_seed_hit_cross_contig_connect/evidence_levels_simple_two_contigs-two_primers_simple.fasta
rm -r test_out_folder

# Test gzipped gff file
call_new_test 'Test gzipped file'
test_exit_status "Magphi -g evidence_levels_single_contig.gff.gz -s two_seeds_w_annotation_between.fasta -o test_out_folder -md 750 -is" 0
test_output_file test_out_folder/master_primer_evidence.csv evidence_level_single_contigs_connect_annotations_include_seeds/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_level_single_contigs_connect_annotations_include_seeds/inter_primer_distance.csv
test_output_file test_out_folder/annotation_primers/evidence_levels_single_contig-annotation_primers.gff evidence_level_single_contigs_connect_annotations_include_seeds/evidence_levels_single_contig-annotation_primers.gff
test_output_file test_out_folder/annotation_primers/evidence_levels_single_contig-annotation_primers.fasta evidence_level_single_contigs_connect_annotations_include_seeds/evidence_levels_single_contig-annotation_primers.fasta
test_output_file test_out_folder/annotation_num_matrix.csv evidence_level_single_contigs_connect_annotations_include_seeds/annotation_num_matrix.csv
rm -r test_out_folder

# Test gzipped fasta file evidence level when only two unique seeds hit with connection but no annotation
call_new_test "Test gzipped fasta file evidence level when only two unique seeds hit with connection but no annotation"
Magphi -g two_contigs_two_primers_single_hit.fasta.gz -s two_primers_simple_match_primers.fasta -o test_out_folder -md 70 > /dev/null 2>&1
test_output_file test_out_folder/master_primer_evidence.csv evidence_levels_simple_hits_cross_contig_connected/master_primer_evidence.csv
test_output_file test_out_folder/inter_primer_distance.csv evidence_levels_simple_hits_cross_contig_connected/inter_primer_distance.csv
test_output_file test_out_folder/two_primers_simple/two_contigs_two_primers_single_hit-two_primers_simple_1_break.fasta evidence_levels_simple_hits_cross_contig_connected/two_contigs_two_primers_single_hit-two_primers_simple_1_break.fasta
test_output_file test_out_folder/two_primers_simple/two_contigs_two_primers_single_hit-two_primers_simple_2_break.fasta evidence_levels_simple_hits_cross_contig_connected/two_contigs_two_primers_single_hit-two_primers_simple_2_break.fasta
rm -r test_out_folder


# Chaws problem.
call_new_test "Test larger mock up of problem reported in emm4"
Magphi -g larger_Magphi_test/*.gff -s larger_test_primers.fa -o test_out_folder -md 15000 -is
test_output_file test_out_folder/annotation_num_matrix.csv larger_Magphi_test/annotation_num_matrix.csv
test_output_file test_out_folder/inter_primer_distance.csv larger_Magphi_test/inter_primer_distance.csv
test_output_file test_out_folder/contig_hit_matrix.csv larger_Magphi_test/contig_hit_matrix.csv
test_output_file test_out_folder/master_primer_evidence.csv larger_Magphi_test/master_primer_evidence.csv
rm -r test_out_folder



# 3. End of testing - check if any errors occurrred
if [ "$num_errors" -gt 0 ]; then
    echo "$test_program failed $num_errors out of $num_tests tests"
    exit 1
else
    echo "$test_program passed all $num_tests successfully"
    exit 0
fi

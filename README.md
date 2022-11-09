[![Test](https://github.com/milnus/Magphi/actions/workflows/test.yml/badge.svg)](https://github.com/milnus/Magphi/actions/workflows/test.yml)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04369/status.svg)](https://doi.org/10.21105/joss.04369)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/magphi/README.html)


# Magphi
**Tool intended to pull out sequences and annotations from genomes between user provided fasta sequence paris (seed sequences) using BLAST**

## When to use
if you have moderately large (100+bp) sequences such as genes that you know flank a particular region of your genome
such as prophage attachment sites and wish to extract the region in one or more genomes. Particularly useful if there is
some level of acceptable variation in the site that will allow a hit with blast but not a direct pattern match.
If you have small (<100 bp) motifs with confidence as to the number of ambiguities and direction of primer this tool
may not be ideal due to the inability for BLAST to deal with short sequences. consider trying [Seqkit](https://bioinf.shenwei.me/seqkit/) instead.

## Installation
### Conda installation
``` conda install -c bioconda magphi ```

### Pip install
Make sure you have the right versions of Samtools and Bedtools installed  
``` pip install Magphi```

#### Dependencies
* Python >= 3.9
* Blast
* Bedtools
* Samtools
* pybedtools


## Help command
```
usage: Magphi [-h] -g .fa/.gff [.fa/.gff ...] -s multi_fasta_file.fa [-is]
              [-md int] [-p] [-S] [-b | -n] [-o path/to/output] [-c int]
              [-l | -q] [--check] [-v]

Welcome to Magphi! This program will extract sequences and possible
annotations within a set of seed sequences given as input.

optional arguments:
  -h, --help            show this help message and exit
  -g .fa/.gff [.fa/.gff ...], --input_genomes .fa/.gff [.fa/.gff ...]
                        Give the fasta or gff3 files. (gff3 files should
                        contain the genome fasta sequence)
  -s multi_fasta_file.fa, --input_seeds multi_fasta_file.fa
                        Give the multi fasta containing the seed sequences to
                        be used for extracting sequnces
  -is, --include_seeds  Argument to include the seeds in the
                        sequence/annotations extracted [default: seeds are
                        removed]
  -md int, --max_seed_distance int
                        The maximum distance with which seeds will be merged.
                        This can often be set a bit higher than an expected
                        size of a region If no maximum distance is wanted then
                        set to 0 [default: 50,000bp]
  -p, --protein_seed    to use tblastn instead of blastn when protein seeds
                        are supplied - useful for hits across diverse genomes
  -S, --stop_orientation
                        Argument to NOT reorient output sequences and
                        annotations by first seed in pair (Only for connected
                        seed not contig breaks) [default: sequences and
                        annotations are oriented]
  -b, --print_breaks    Argument to print outputs when seeds are next to
                        contig breaks [default: sequences are not printed]
  -n, --no_sequences    Argument to not print outputs related to annotations
                        or sequences found between seeds [default: sequences
                        are printed]
  -o path/to/output, --output_folder path/to/output
                        Give path to output folder [default: current folder]
  -c int, --cpu int     Give max number of CPUs [default: 1]
  -l, --log             Record program progress in for debugging purpose
  -q, --quiet           Only print warnings
  --check               Check dependencies for Magphi and exit
  -v, --version         show program's version number and exit

  ```

## Constructing seed sequence file
Seed sequences are what define a region in which the user is interested. These sequences can be core genes, or known stretches of DNA that surround a specific region in genomes. It is recommended that seed sequences are no less than 100 bp, as BLAST does not work well for smaller pieces of DNA.

### Rules to follow for seed sequences
* Seed sequences should be no less than 100 bp
* Names should be unique for each pair of seed sequences
* Each name of a seed sequence in a pair should be appended with a unique tag (i.e. ```_1``` or ```_2```) to make them unique within the pair.
* If a seed sequences is used twice it should be duplicated in the fasta file and given two unique names

### Example seed sequence file:
```
>phageA_1
CTGGTCGGTTGGTATATAGTGTTCGCACTG.....
>phageA_2
CACTTGTAACACAACCAGACCCCCCCGGAA.....
>spec_phage1
CACTTGTAACACAACCAGACCCCCCCGGAA.....
>spec_phage2
TGCAGTTCCCCATGGGGTGTAGGTGTAACG.....
```

## Genomes
Genomes should be in either Fasta format or GFF3 format with the genome appended in the end of the file, with the ##FASTA line separating the annotations from the appended fasta genome.
Only a single file format (either Fasta or GFF) is allowed in a single run. Files can be Gzipped, but a mix of Gzipped and non-Gzipped files are not allowed.

### Rules to follow for input genomes
* Fasta or GFF3 format with appended genome
* Exclusively Fasta or GFF files for a single run, no mixing
* Genomes can be Gzipped, but only Gzipped or not for a single run.

## Results
For each seed sequence pair given as input a subfolder in the output folder will be created. In a seed sequence subfolder extracted genomic sequences and annotations can be found.

Additional outputs are tables summarising findings from the Magphi run. Besides ```seed_pairing.tsv``` the output tables all contain the name of genomes given as the first column. The subsequent column headers are the seed sequences identified for the run.

### Summary tables
* ```seed_pairing.tsv``` - Gives the common name identified for a seed sequence pair and the two seed sequences found in the pair. Can be used to check if Magphi identified seed sequence pairs as intended by the user.
* ```contig_hit_matrix.csv``` - Indicates the number of locations a pair of seed sequences were identified to map in a given genome by BLAST. Can be used to evaluate seed specificity.
* ```inter_seed_distance.csv ``` - Genomic distance in basepairs between merged seed sequences in a pair. Additional columns will be added for each seed sequence that is found to be next to a sequence break. - Can be used to evaluate insertions of genetic material between seed sequences or varying distance between genomic features.
* ```annotation_num_matrix.csv``` - The number of annotated features extracted from a GFF file, if any can be found between merged seed sequences. This can be used to evaluate changes in gene/genetic feature content between seed sequences. This is useful for searching insertion sites of mobile genetic elements.
* ```master_seed_evidence.csv``` - Overview of how each seed sequence pair was found to map, merge, and extract features for a specific genome. This is easily the most useful output for evaluating seed sequence quality, and large scale impression of output quality.

### Outputs to expect
The standard outputs listed above are to be expected for all Mapghi runs. 
Additional outputs related to extracted fasta sequences and gff annotations can be expected in some instances. When running with:
* Default parameters (no ```-b``` or ```-n``` arguments to Magphi), Fasta outputs can be expected for evidence levels 5B and 5C. A Gff output can be expected with evidence level 5C.
* When ```-b``` is given to Magphi, Fasta outputs can be expected for evidence levels mentioned above, and 4B and 4C. Gff outputs can be expected for 4C as well.
* No output is expected when ```-n``` is given as a parameter.
* As a default output fasta and Gff files are oriented so that the first seed in a pair (see ```seed_pairing.tsv``` file) is at the beginning of the sequence on the positive strand

```-n``` and ```-b``` are mutually exclusive as they both alter the expected output and therefore Magphi does not allow both to be given in a single command.

### Evidence level overview
0. One or no seed sequence hit the genome
1. Multiple seed sequences hit with no overlap
2. Multiple seed sequences hit with multiple connections
3. Seed(s) deleted due to overlap or placed at end of contig and seeds are excluded
4. 
   A. Two Seed sequences hit on seperate contigs - No connection

   B. Two Seed sequences hit on seperate contigs - with connection no annotations

   C.Two Seed sequences hit on seperate contigs - with connection and annotations
5. 
   A. Two seed sequences hit on same contig - not connected

   B. Two seed sequences hit on same contig - connected no annotations
   
   C. Two seed sequences hit on same contig - connected with annotations

## Improvements
Depending on the evidence level some changes can be made to improve a Magphi run. Below are some points to consider:
* Evidence level of 0 can indicate a bad selection of seed sequences or poor quality of assembly.
* Too small maximum distance or poor seed sequences can result in an evidence level of 1
* Too large a maximum distance may result in draft genomes having an evidence level of 2, however decreasing the max distance may result in false region being extracted
* When 3 is given as evidence level try ```-ip``` to see the resulting evidence level, to examine if seed sequences overlap or fall on edge of contig.
* Evidence level = 4 can give output if ```-b``` is given as an argument to Mapghi.
* Dividing seed sequence pairs with similar max distance into separate runs is a good trick to maximise likelihood of good extractions.
* For large datasets consider using ```-n``` to save memory, if your analysis is not interested in the sequence itself.


## Community contributions via Github
All questions via the [Issues tab](https://github.com/milnus/Magphi/issues) are answered as quickly as possible!
 - The more information you can give about problems and/or enquiries the easier it is to help you quickly. Aim to narrow down the problem and give a way to reproduce the problem, if any.    

Any suggestions to improvements are welcome!
 - Feel free to suggest new features or ways of doing things better. Things will be implemented as good as possible

Any contributions are appreciated! 
 - If you have the time and will feel free to send a pull request. Unit- and functional tests are not required. Try to comment the code in a way that allows for easier understanding. Use doc-stings for each function to explain their in- and output.

## For more info
See Wiki tab for more info on the workings, inputs, and outputs of Magphi

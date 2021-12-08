[![Test](https://github.com/milnus/Magphi/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/milnus/Magphi/actions/workflows/test.yml)

# Magphi
Tool intended to pull out sequences and annotations between seed sequences specified by the user.

## When to use
# TODO!

## Installation
### Conda installation
``` # Comming soon ```

### Pip install
Make sure you have the right versions of Samtools and Bedtools installed

``` pip install Magphi```

### Setting up conda environment using Conda 
Solution until conda recipie has been created  
```$ conda create -n Magphi -c bioconda samtools=1.13 blast bedtools python=3.9```  
```$ conda activate Magphi```  
```$ pip install Magphi```  
To test the insatllation of dependencies run:  
```$ Magphi --check```

## Constructing seed seqeunce file
Seed sequences are what define a region in which the user is interested. These seqeunces can be core genes, or known stretches of DNA that surround a specific region in genomes. It is recommended that seed seqeunces are no less than 100 bp, as BLAST does not work well for smaller pieces of DNA.

### Rules to follow for seed sequnces
* Seed seqeunces should be no less than 100 bp
* Names should be unique for each pair of seed sequnces
* Each name of a seed seqeunce in a pair should be appended with a unique tag (i.e. ```_1``` or ```_2```) to make them unique within the pair.
* If a seed seqeunces is used twice it should be duplicated in the fasta file and given two unique names

### Examlpe seed seqeunce file:
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
Genomes should be in either Fasta format or GFF3 format with the genome appended in the end of the file, with the ##FASTA line seperating the annotations from the appended fasta genome.
Only a single file format (either Fasta or GFF) is allowed in a single run. Files can be Gzipped, but a mix of Gzipped and non-Gzipped files are not allowed.

### Rules to follow for input genomes
* Fasta or GFF3 format with appended genome
* Exclusively Fasta or GFF files for a single run, no mixing
* Genomes can be Gzipped, but only Gzipped or not for a single run.

## Results
For each seed seqeunce pair given as input a subfolder in the output folder will be created. In a seed seqeunce subfolder extracted genetic seqeunces and annotations can be found.

Additional outputs are tables summarising findings from the Magphi run. Besides ```seed_pairing.tsv``` the output tables all contain the name of genomes given as the first column. The subsequent column headers are the seed sequences identified for the run.

### Summary tables
* ```seed_pairing.tsv``` - Gives the common name identified for a seed seqeunce pair and the two seed seqeunces found in the pair. Can be used to check if Magphi identified seed seqeunce pairs as intended by the user.
* ```contig_hit_matrix.csv``` - Indicates the number of locations a pair of seed seqeunces were identifed to map in a given genome by BLAST. Can be used to evaluate seed specificity.
* ```inter_seed_distance.csv ``` - Genomic distance in basepairs between merged seed seqeunces in a pair. Additional columns will be added for each seed seqeunce that is found to be next to a seqeunce break. - Can be used to evaluate insertions of genetic material between seed seqeuences or varying distance between genomic features.
* ```annotation_num_matrix.csv``` - The number of annotated features extracted from a GFF file, if any can be found between merged seed seqeunces. This can be used to evaluate changes in gene/genetic feature content between seed seqeunces. This is usefull for searching insertion sites of mobile genetic elements.
* ```master_seed_evidence.csv``` - Overview of how each seed seqeunce pair was found to map, merge, and extract features for a specific genome. This is easily the most usefull output for evaluating seed seqeunce quality, and large scale impression of output quality.

### Outputs to expect
The standard outputs listed above are to be expected for all Mapghi runs. 
Additional outputs related to extracted fasta seqeunces and gff annotations can be expected in some instances. When running with:
* Default parameters (no ```-b``` or ```-n``` arguments to Magphi), Fasta outputs can be expected for evidence levels 5B and 5C. A Gff output can be expected with evidence level 5C.
* When ```-b``` is given to Magphi Fasta outputs can be expected for evidence levels mentioned above, and 4B and 4C. Gff outputs can be expected for 4C as well.
* No output is expected when ```-n``` is given as a parameter.

```-n``` and ```-b``` are mutually exclusive as they both alter the expected output and therefor Magphi does not allow both to be given in a single command.

### Evidence level overview
0. One or no seed seqeunce hit the genome
1. Multiple seed seqeunces hit with no overlap
2. Multiple seed seqeunces hit with multiple connections
3. Seed(s) deleted due to overlap or placed at end of contig and seeds are excluded
4. 
   A. Two Seed sequences hit on seperate contigs - No connection

   B. Two Seed sequences hit on seperate contigs - with connection no annotations

   C.Two Seed sequences hit on seperate contigs - with connection and annotations
5. 
   A. Two seed seqeunces hit on same contig - not connected

   B. Two seed seqeunces hit on same contig - connected no annotations
   
   C. Two seed seqeunces hit on same contig - connected with annotations

## Improvements
Depending on the evidence level some changes can be made to imporive a Magphi run. Below are some points to consider:
* Evidence level of 0 can indicate a bad selection of seed seqeunces or poor quality of assembly.
* Too small maxium distance or poor seed seqeunces can result in an evidence level of 1
* Too large a maximum distance may result in draft genomes having an evidence level of 2, however decreasing the max distance may result in false region being extracted
* When 3 is given as evidence level try ```-ip``` to see the resulting evidence level, to examine if seed sequences overlap or fall on edge of contig.
* Evidence level = 4 can give output if ```-b``` is given as an argument to Mapghi.
* Dividing seed sequence pairs with similar max distance into seperate runs is a good trick to maximise likelihood of good extractions.
* For large datasets consider using ```-n``` to save memory, if your analysis is not interested in the seqeunce itself.

## For more info
See Wiki tab for more info on the workings, inputs, and outputs of Magphi


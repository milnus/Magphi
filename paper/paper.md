---
title: 'Magphi: **Automated** sequence and feature extraction from genomes using seed sequences'
tags:  
  - Python
  - microbiology
  - genomics
  - genome
authors:
  - name: Magnus G. Jespersen 
    orcid: 0000-0001-9751-9877
    affiliation: 1
  - name: Andrew Hayes
    orcid: 0000-0001-8038-1656
    affiliation: 1
  - name: Mark R. Davies
    orcid: 0000-0001-6141-5179
    affiliation: 1
affiliations:
  - name: Department of Microbiology and Immunology, University of Melbourne at the Peter Doherty Institute for Infection and Immunity, Melbourne, VIC, Australia
    index: 1
date: 13 October 2021  
bibliography: paper.bib
---

# Summary
With the decrease in cost of sequencing genomes from biological organisms the burden of work has been moved from the generation of genome sequences to the analytics. Researchers working with genomes originating from microorganisms often work with multiple genomes in a single analysis. The number of genomes in datasets can lead to challenges when it comes to extracting specific regions of interest from multiple genomes. Manual extraction of regions becomes impractical and time consuming when datasets exceed 10-20 genomes. The complexity of this task increases when working within complex regions of genomes that may not assemble into a single contiguous region using some existing technologies such as short read-based sequencing technologies. Therefore, **automation** is required as datasets of microbial genomes routinely consist of tens or hundreds of genomes. Here we present Magphi, a BLAST [@mount:2007] based  contig-aware genome extraction tool utilising seed sequences to identify and extract regions of interest.

# Statement of need
Magphi aims to free researchers time from tedious tasks and allow them to apply domain knowledge and focus on analysing genomic regions of interest. Packages such as Seqkit [@shen:2016] allows for extraction and manipulation of Fasta and Fastq files; However, such tools are not aimed at working with Gene Feature Format (GFF) files, a common file type in bioinformatics, or when regions of interest may span across contigs. Handling of Gene Feature Format (GFF) files are often necessary when researchers examine genomes annotated genes and other genetic features, as these are not included in Fasta formatted files.  

Magphi is a command-line tool written in Python 3. It uses the Basic Local Alignment Search Tool (BLAST) [@mount:2007], BEDtools [Quinlan:2010] and implements logic to identify possible connections between given seed-sequences to return the optimal solution in terms of genetic sequence and possible annotations between a set of seed-sequences. Magphi can handle fasta files or GFF version 3 files with included genomes, as given by the microbial annotations tool Prokka [@seemann:2014]. Magphi is contig-aware, and will return a file containing confidence level for each pair of seed-sequences and genomes, providing the researcher with fine grained feed-back on their run. The file containing confidence levels can be imported into Phandango [@Hadfield:2017], along with a phylogenetic tree of genomes for quick and visual inference of patterns of potential problems.

Magphi is scalable and can as input take multiple genomes and multiple pairs of seed-sequences. Outputs are divided by the input seed-sequences for easier file management

# Acknowledgements
MGJ is supported by The Melbourne Research Scholarship from The University of Melbourne

# References
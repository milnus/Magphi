---
title: 'Magphi: Sequence extraction tool from FASTA and GFF3 files using seed pairs'
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
date: 16 December 2021  
bibliography: paper.bib
---

# Summary
Researchers working with genomes originating from microorganisms often work with multiple genomes in a single analysis. The number of genomes in datasets can pose challenges when it comes to extracting specific regions of interest from multiple genomes. Manual extraction of regions becomes impractical and time consuming when datasets exceed 10-20 genomes. The complexity of this task increases when working within complex regions of genomes that may not assemble into a single contiguous sequence using some existing technologies such as short read-based sequencing technologies. Therefore, automation is required as datasets of microbial genomes routinely consist of tens or hundreds of genomes. Here we present Magphi, a BLAST [@mount:2007] based  contig aware genome extraction tool utilising seed sequences to identify and extract regions of interest.

# Statement of need
Magphi extracts genomic regions of interest from FASTA and Gene Feature Format 3 (GFF3) files, both being common file types in bioinformatics. Packages such as Seqkit [@shen:2016] allow for extraction and manipulation of FASTA and FASTQ files; However, such tools do not work with GFF3, or when regions of interest may span across contigs. Handling of GFF3 files are often necessary when researchers examine annotated genomes, as these are not included in FASTA formatted files.  

Magphi is a command-line tool written in Python 3. It uses the Basic Local Alignment Search Tool (BLAST) [@mount:2007], BEDtools [@Quinlan:2010] and implements logic to identify possible connections between given seed sequences to return the optimal solution in terms of genetic sequence and possible annotations between a set of seed sequences. Magphi can handle FASTA or GFF 3 files with included genomes, as given by the microbial annotation tool Prokka [@seemann:2014]. Magphi is contig aware, and will return a file containing the confidence level for each pair of seed sequences and genomes, providing the researcher with feedback on their run. The file containing confidence levels, distances between seed sequences, and number of annotations can be imported into Phandango [@Hadfield:2017], along with a phylogenetic tree of genomes for quick and visual inference of patterns or potential problems. Magphi also produces an output folder for each seed sequence pair, containing FASTA and GFF3 files when possible.   

Magphi is scalable and can take multiple genomes and pairs of seed sequences. Outputs are divided by the input seed sequences for easier file management

# Acknowledgements
MGJ is supported by The Melbourne Research Scholarship from The University of Melbourne.

# References
---
title: 'Magphi: Automated sequence and feature extraction from genomes using seed sequences'
tags:  
- Python
- microbiology
- genomics
- genome
- sequence extraction
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
- Department of Microbiology and Immunology, University of Melbourne at the Peter Doherty Institute for Infection and Immunity, Melbourne, VIC, Australia
  index: 1
date: 13 October 2021  
bibliography: paper.bib
---

# Summary
With the decrease in cost of sequncing genomes from biological organisms the burden of work has been moved from the generation of genome sequences to the analysis. Researchers working with genomes originating from micro organism often work with multiple genomes in a single analysis. The number of genomes in datasets can lead to challenges when it comes to extracting specific regions of interest of genomes from multiple genomes. Manual extraction of regions becomes impractical and timeconsuming when datasets exceed 10-20 genomes. The complexity of this task increases with some regions of genomes not being well assembled when reconstructing a genome from some current technologies (Illumina short reads sequencing). Therefore, automation is required as datasets of microbial genomes routinely consist of tenths or hundreds of genomes, with more and more having thousands of genomes [ref??]. Here we Magphi, a BLAST [@mount:2007] based  contig-aware genome extraction tool utilising seed sequences to identify regions of interest.

# Statement of need
Magphi aims to free researchers time from tedious tasks and allow them to apply domain knowledge and focus on analysing results. New packages within genome research address the same issue, one being Seqkit [@shen:2016]. Seqkit is a manipulation tool aimed at standard file formats in bioinformatics, fasta and fastq. Seqkit allows for the extraction of genome sequences from fasta files, using multiple functions. However, tools for working with Gene Feature Format (GFF) files are required, when researchers handle genomes annotated with genes and other genetic features.  

Magphi build on the Basic Local Alignment Search Tool (BLAST) [@mount:2007] and implements logic to identify possible connections between given seed-seqeunces to return the optimal solution in terms of genetic sequence and possible annotations between a set of seed-seqeunces. Magphi can handle fasta files or GFF version 3 files with included genomes, as given by the microbial annotations tool Prokka [@seemann:2014]. Magphi is contig-aware, and will return a file containing confidence level for each pair of seed-seqeunces and genomes, providing the researcher with fine grained feed-back on their run. The file conteining confidence levels can be imported into Phandango [@Hadfield:2017], along with a phylogenetic tree of genomes for quick and visual inference of patterns of potential problems.

Magphi is scalable and can as input take multiple genomes and multiple pairs of seed-seqeunces. Outputs are divided by the input seed-seqeune for easier file management

# Acknowledgements
MGJ is supported by The Melbourne Research Scholarship from The University of Melbourne

# References
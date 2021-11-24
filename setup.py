#!/usr/bin/env python

from setuptools import setup

LONG_DESCRIPTION = \
'''The program extracts regions of interest from Fasta or Genome Feature Format (GFF) genomes.
This is done given a set of seed sequences given as nucleotide strings in a multi-line fasta file.  
The program can output fasta and GFF outputs or regions, and will giv multiple outputs around regions and their evidence.

The program takes in multiple fasta or GFF files at a single time and can also take multiple seed sequence pairs.'''


setup(
    name='Magphi',
    version='0.1.3',
    author='Magnus Ganer Jespersen',
    author_email='magnus.ganer.j@gmail.com',
    packages=['Magphi'],
    package_dir={'Magphi': 'Magphi'},
    entry_points={
        'console_scripts': ['Magphi = Magphi.__main__:main']
    },
    url='https://github.com/milnus/Magphi',
    license='MIT license',
    description=('A bioinformatics tool allowing for examnination and extraction of genomic features using seed sequences.'),
    long_description=LONG_DESCRIPTION,
    install_requires=['biopython==1.79',
                      'pybedtools'],
    keywords=['Genomic', 'extraction', 'bacteria', 'prokaryotes', 'bioinformatics'],
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 4 - Beta']
)

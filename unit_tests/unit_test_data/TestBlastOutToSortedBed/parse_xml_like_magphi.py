from Bio import SearchIO
import pybedtools as bedtools
blast_xml_output = open('Mock_blast_out.xml', "r") # open blast xlm file as object
print( blast_xml_output.read())
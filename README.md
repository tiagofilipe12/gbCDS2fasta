# gbCDS2fasta.py

This script retrieves CDS data from .gb files to a fasta in which each CDS is a sequence entry. Currently it supports .gb files as input files. Fasta files can also be provided for sequence retrieving if it perfectly matches the references in .gb files.

**'-in'**,**'--input'**, dest='inputfile', required=True, help='Provide the input directory with .gb (GenBank) files to parse.'

**'-f'**,**'--fasta'**, dest='fastahelper', nargs='+', help='If, for some reason, your .gb file does not have the nucleotide sequence data you need, you can use this option. It allows to provide a fasta with the same accession numbers as the .gb file that you are trying to converto into single CDS fasta sequences.'

**'-out'**,**'--output'**, dest='outputfile', required=True, help='Provide the output file name.'

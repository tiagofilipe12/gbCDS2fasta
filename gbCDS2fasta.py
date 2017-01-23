#!/usr/bin/env python

## Last update: 23/1/2017
## Author: T.F. Jesus
## This script retrieves CDS data from .gb files to a fasta in which each CDS is a sequence entry
## Untested script and lacks some functions

import argparse
from Bio import SeqIO

def header_fix(input_header):
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/"]
	for char in problematic_characters:
		input_header=input_header.replace(char, '_')
	return input_header

def genbank_file_reader(input_dir, output_file):	
	for infile in os.listdir(input_dir):
		if input_dir.endswith(".gb"):
			cds_gene_position = {}	
			record = SeqIO.read(infile, "genbank")
			gb_description=gbkdata.description[:-1]
			seq_description=header_fix(gb_description)
			sequence = str(gbkdata.seq)			
			seq_id = record.id
			tag_name = "ref_" + seq_id + "_" + seq_description + "_"  
			for feat in record.features:
				strand = feat.location.strand
				seq_type = feat.type
				if seq_type == "CDS":
					location = feat.location ## did not use .start and .end because it does not retrieve properly when DNA sequence is circular
					if "gene" in feat.qualifiers.keys():
						tag_name =+ "CDS_" + feat.qualifiers.values()[0]
						cds_gene_position[tag_name] = location
			cds2fasta(output_file, cds_gene_position, sequence)	

def cds2fasta(output_file, cds_gene_position, sequence):
	out_handle=open(output_file,'w')
	for k in cds_gene_position.keys():
		out_handle.write(">" + cds_gene_position.keys() + "\n")
		out_handle.write(sequence)

def main():
	parser.add_argument('-in','--input', dest='inputfile', required=True, help='Provide the input fasta files to parse')
	parser.add_argument('-out','--output', dest='outputfile', required=True, help='Provide the output directory')
	args = parser.parse_args()
	list_of_gb = args.inputfile.split(" ")
	for gb in list_of_gb:
		genbank_file_reader(gb, args.outputfile)

if __name__ == "__main__":
	main()
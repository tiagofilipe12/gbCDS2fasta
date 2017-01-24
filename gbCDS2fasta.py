#!/usr/bin/env python

## Last update: 24/1/2017
## Author: T.F. Jesus
## This script retrieves CDS data from .gb files to a fasta in which each CDS is a sequence entry
## Untested script and lacks some functions

import argparse
import os
from Bio import SeqIO

def header_fix(input_header):
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/"]
	for char in problematic_characters:
		input_header=input_header.replace(char, '_')
	return input_header

def reverse_complement(sequence):
    rc_dictionary = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([rc_dictionary[base] for base in reversed(sequence)])

def genbank_file_reader(input_dir, output_file):	
	out_handle=open(output_file,'w')
	for infile in os.listdir(input_dir):
		if infile.endswith(".gb"):
			cds_gene_position = {}	
			record = SeqIO.read(os.path.join(input_dir,	infile), "genbank")
			gb_description=record.description[:-1]
			seq_description=header_fix(gb_description)
			sequence = str(record.seq)			
			seq_id = record.id
			tag_name = "ref_" + seq_id + "_" + seq_description + "_"  
			for feat in record.features:
#				strand = feat.location.strand
				seq_type = feat.type
				if seq_type == "CDS":
					cds_tag_name = ""
					location = str(feat.location) ## did not use .start and .end because it does not retrieve properly when DNA sequence is circular
					for key in feat.qualifiers.keys():
						cds_tag_name = tag_name + "CDS_" + location
					cds_gene_position[cds_tag_name] = location
			cds2fasta(out_handle, cds_gene_position, sequence)	

def cds2fasta(out_handle, cds_gene_position, sequence):
	counter_start = 0
	counter_end = 0
	for k in cds_gene_position.keys():
		print k		
		start_location = cds_gene_position[k].split(":")[0].lstrip("[")
		end_location = cds_gene_position[k].split(":")[1].split("]")[0]
		complement = cds_gene_position[k].split(":")[1].split("]")[1] ## (+) or (-)
		if ">" or "<" in start_location:
			counter_start =+1
			start_location = start_location.lstrip("<").lstrip(">")
		if ">" or "<" in end_location:
			counter_end =+ 1
			end_location = end_location.lstrip("<").lstrip(">")
		out_handle.write(">" + k + "\n")
		if "(-)" in complement:
			out_handle.write(reverse_complement(sequence[int(start_location):int(end_location)]) + "\n")
		elif "(+)" in complement: 
			out_handle.write(sequence[int(start_location):int(end_location)] + "\n")
	return counter_start, counter_end

def main():
	parser = argparse.ArgumentParser(description="Converts all CDS in a .gb file into a fasta with each CDS region as an individual sequence")
	parser.add_argument('-in','--input', dest='inputfile', required=True, help='Provide the input directory with .gb (GenBank) files to parse')
	parser.add_argument('-out','--output', dest='outputfile', required=True, help='Provide the output file name')
	args = parser.parse_args()
	list_of_gb = args.inputfile.split(" ")
	for gb in list_of_gb:
		genbank_file_reader(gb, args.outputfile)

if __name__ == "__main__":
	main()
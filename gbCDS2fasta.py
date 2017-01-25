#!/usr/bin/env python

## Last update: 24/1/2017
## Author: T.F. Jesus
## This script retrieves CDS data from .gb files to a fasta in which each CDS is a sequence entry


import argparse
import os
from Bio import SeqIO

def header_fix(input_header):
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/","[","]",":","{","}"]
	for char in problematic_characters:
		input_header=input_header.replace(char, '_')
	return input_header

def reverse_complement(sequence):
    rc_dictionary = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'X':'X', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V'}
    return "".join([rc_dictionary[base] for base in reversed(sequence)])

def sequencesfromfasta(fasta_files):
	fastasequence_dict = {}
	for fasta in fasta_files:
		if any (x in fasta for x in [".fas",".fasta",".fna",".fsa", ".fa"]):
			fasta_handle=open(fasta,'r')			
			x = 0
			for line in fasta_handle:
				header_split = line.split("|")
				if len(line) > 0: 
					line = line.splitlines()[0] 
				if x==0 and line.startswith(">"):
					accession = header_split[3].strip()
					sequence_list = []
					x+=1
				elif x==0 and not line.startswith(">"):
					print "Is this a fasta file? " + fasta
					raise SystemExit
				elif x >=1 and line.startswith(">"):
					if accession not in fastasequence_dict.keys():
						fastasequence_dict[accession]="".join(sequence_list)
					else:
						print "Something is terribly wrong! More than one hit found for: " + accession
					accession = header_split[3].strip()
					sequence_list=[]
					x=+1
				else:
					sequence_list.append(line)
			if accession not in fastasequence_dict.keys():
				fastasequence_dict[accession]="".join(sequence_list)
			else:
				print "Something is terribly wrong! More than one hit found for: " + accession
			fasta_handle.close()
			print 
		else:
			print "Error: You should provide a fasta file to help constructing CDS cat_sequences. File extension not recognized."
	print str(len(fastasequence_dict.keys())) + " fasta sequences processed."
	return fastasequence_dict    

def genbank_file_reader(input_dir, output_file, fastasequence_dict):	
	out_handle=open(output_file,'w')
	for infile in os.listdir(input_dir):
		if infile.endswith(".gb"):
			cds_gene_position = {}	
			try:
				record = SeqIO.read(os.path.join(input_dir,	infile), "genbank")
			except:
				print infile + " does not have any records on it."
			gb_description=record.description[:-1]
			seq_description=header_fix(gb_description)
			sequence = str(record.seq)
			seq_id = record.id.strip()			
			tag_name = "ref_" + seq_id + "_" + seq_description + "_"  
			for feat in record.features:
				seq_type = feat.type
				if seq_type == "CDS":
					cds_tag_name = ""
					location = str(feat.location) ## did not use .start and .end because it does not retrieve properly when DNA sequence is circular
					for key in feat.qualifiers.keys():
						cds_tag_name = tag_name + "CDS_" + location
					cds_gene_position[cds_tag_name] = location
			counter_start, counter_end = cds2fasta(out_handle, cds_gene_position, sequence, fastasequence_dict, seq_id)
	return counter_start, counter_end

def cds2fasta(out_handle, cds_gene_position, sequence, fastasequence_dict, seq_id):
	counter_start = 0
	counter_end = 0
	fasta_instance = 0
	for k in cds_gene_position.keys():
		cat_sequences = []
		initial_trim = cds_gene_position[k].lstrip("join{").rstrip("}")
		for entry in initial_trim.split(","):
			start_location = entry.split(":")[0].lstrip(" ").lstrip("[")
			end_location = entry.split(":")[1].split("]")[0]
			complement = entry.split(":")[1].split("]")[1] ## (+) or (-)
			if ">" or "<" in start_location:
				counter_start =+1
				start_location = start_location.lstrip("<").lstrip(">")
			if ">" or "<" in end_location:
				counter_end =+ 1
				end_location = end_location.lstrip("<").lstrip(">")
			if sequence.startswith("N"):
				sequence = fastasequence_dict[seq_id]
			else:
				pass
			if "(-)" in complement:
				cat_sequences.append(reverse_complement(sequence[int(start_location):int(end_location)]))
			elif "(+)" in complement:
				cat_sequences.append(sequence[int(start_location):int(end_location)])
		out_handle.write(">" + header_fix(k) + "\n")
		out_handle.write("".join(cat_sequences) + "\n")
	return counter_start, counter_end

def main():
	parser = argparse.ArgumentParser(description="Converts all CDS in a .gb file into a fasta with each CDS region as an individual sequence")
	parser.add_argument('-in','--input', dest='inputdir', nargs='+', required=True, help='Provide the input directory with .gb (GenBank) files to parse')
	parser.add_argument('-f','--fasta', dest='fastahelper', nargs='+', help='If, for some reasion, your .gb file does not have the sequence data you need, you can use this option that allow to provide a fasta with the same accession numbers of the entries on the.gb file you are trying to convert into single CDS fastas')
	parser.add_argument('-out','--output', dest='outputfile', required=True, help='Provide the output file name')
	args = parser.parse_args()
	list_of_gb = args.inputdir
	list_of_fasta = args.fastahelper
	if list_of_fasta:
		fastasequence_dict = sequencesfromfasta(list_of_fasta)
	else:
		fastasequence_dict = 0
	for gb in list_of_gb:
		counter_start, counter_end = genbank_file_reader(gb, args.outputfile, fastasequence_dict)
	print str(counter_start) + " sequences with starting fuzzy positions. Please refer to http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc36"
	print str(counter_end) + " sequences with ending fuzzy positions. Please refer to http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc36"
	print "Note: This script does not handle fuzzy locations, rather it takes the location provided by NCBI record for this position."

if __name__ == "__main__":
	main()

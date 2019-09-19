#!/usr/bin/env python
##########################################################################
# Function: build reference index and annotation info file
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Version: Version 1.0	07/14/2018
##########################################################################
import sys
import os
import re
import argparse

#==============================================================
# index with bowtie2
#==============================================================
def bowtie2_index (ref_seq, base_name):
	print(">>Running botwie2 to index reference...")
	# Index reference sequences
	print("bowtie2-build " + ref_seq + " " + base_name + "\n")
	os.system("bowtie2-build " + ref_seq + " " + base_name)
# function bowtie2_index


#==============================================================
# Build SAF (simplified annotation format) info
#==============================================================
def extract_SAF_info_gene (seq_file, bed_out):	# use NR gene category as the reference 
	open_in = open(seq_file, "r")
	outs = {}
	myid = ""
	for line in open_in.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):	# new gene
			mym = re.search(">([\S]+)", line)
			myid = mym.group(1)
			outs[myid] = 0
			continue
		else:
			outs[myid] = outs[myid] + len(line)
	# foreach line
	open_in.close()
	
	open_out = open(bed_out, "w")
	open_out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
	for myid in sorted(outs.keys()):
		open_out.write(myid + "\t" + myid + "\t1\t" + str(outs[myid]) + "\t+\n")
	# foreach gene
# extract_SAF_info_gene


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-r', help='the reference file', required=True)
	parser.add_argument('-b', help='the reference index base name', required=True)
	parser.add_argument('-t', help='the type of reference, e.g. gene | contig', required=True)
	parser.add_argument('-o', help='the output file for simplified annotation file', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start gene_abundance_indexRef.py -r " + values.r + " ####\n")
	
	### Index ###
	sys.stderr.write("Bowtie2 index......starting\n")
	bowtie2_index (values.r, values.b) 
	sys.stderr.write("Bowtie2 index......done\n")
	
	### Getting SAF ###
	sys.stderr.write("\nExtract SAF info......starting\n")
	if values.t == "gene":
		extract_SAF_info_gene (values.r, values.o)
	sys.stderr.write("\nExtract SAF info......done\n")

	sys.stderr.write("### Finish gene_abundance_indexRef.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()
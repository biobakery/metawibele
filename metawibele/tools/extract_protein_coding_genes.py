#!/usr/bin/env python
##########################################################################
# Function: Extract protein coding genes from total gene sets
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Version: Version 1.0	07/10/2018
##########################################################################
import sys
import os
import re
import argparse

#==============================================================
# collect protein coding IDs
#==============================================================
def collect_seq_id (aa_seq):	#combined_peptides.sorted.faa  
	seqid = {}
	open_file = open(aa_seq, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue                      
		if not re.search("^>", line):
			continue
		mym = re.match(">([\S]+)", line)
		myid = mym.group(1)
		seqid[myid] = myid
	# foeach line
	open_file.close()
	return  seqid
# function collect_seq_id

#==============================================================
# output sequence info
#==============================================================
def output_info (fna_seq, seqid, outfile):
	open_in = open (fna_seq, "r")
	open_out = open(outfile, "w")
	flag = 0
	for line in open_in:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.match(">([\S]+)", line)
			myid = mym.group(1)
			if myid in seqid:
				flag = 1
			else:
				flag = 0
		if flag == 1:
			open_out.write(line + "\n")
	# foreach line
	open_in.close()
	open_out.close()
	#foreach line
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-g', help='input total gene nucleotide sequences file', required=True)
	parser.add_argument('-p', help='input protein coding AA sequences file', required=True)
	parser.add_argument('-o', help='output protein coding nucleotide sequence file', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start extract_protein_coding_genes.py -g " + values.g + " ####\n")
	
	### collect sequence id info ###
	sys.stderr.write("Get seq id info ......starting\n")
	seqid = collect_seq_id (values.p)
	sys.stderr.write("Get seq id info ......done\n")
	
	### Output protein coding gene nuclotide sequence
	sys.stderr.write("\nOutput sequence info ......starting\n")
	output_info (values.g, seqid, values.o)
	sys.stderr.write("Output sequence info ......done\n")

	sys.stderr.write("### Finish extract_protein_coding_genes.py ####\n\n\n")

# end: main

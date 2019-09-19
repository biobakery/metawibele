#!/usr/bin/env python
##########################################################################
# Function: Extact representative amino acid sequences from the clusters based on nucleotide sequences
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Version: Version 1.1	07/01/2018
##########################################################################
import sys
import os
import re
import argparse

#==============================================================
# collect sequences
#==============================================================
def collect_seq_id (nr_seq):	# combined_genes_meta_PC.sorted.clust.rep.fna
	seqid = {}
	open_file = open(nr_seq, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue                      
		if re.search("^>", line):
			mym = re.search(">([\S]+)", line)
			myid = mym.group(1)
			seqid[myid] = ""
	# foeach sample
	open_file.close()
	return  seqid
# function collect_seq_id

#==============================================================
# output sequence info
#==============================================================
def output_info (AA_seq, seqid, outfile):
	open_in = open (AA_seq, "r")
	open_out = open(outfile, "w")
	flag = 0
	flags = {}
	for line in open_in.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([\S]+)", line)
			myid = mym.group(1)
			if myid in seqid:
				flag = 1
				flags[myid] = ""
			else:
				flag = 0
		if flag == 1:
			open_out.write(line + "\n")
	# foreach line
	open_in.close()
	open_out.close()

	# check sequences
	for myid in seqid.keys():
		if not myid in flags:
			print("No AA sequences for gene_id\t" + myid)
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-n', help='input nr nucleotide file', required=True)
	parser.add_argument('-f', help='input original AA file', required=True)
	parser.add_argument('-o', help='output non-redundance AA sequence file', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start extract_non_redundance_AA_seq.py -n " + values.n + " ####\n")
	
	### collect cluster id info ###
	sys.stderr.write("Get clustering info ......starting\n")
	seqid = collect_seq_id (values.n)
	sys.stderr.write("Get clustering info ......done\n")
	
	### Output non-redundance AA sequence
	sys.stderr.write("\nOutput sequence info ......starting\n")
	output_info (values.f, seqid, values.o)
	sys.stderr.write("Output sequence info ......done\n")

	sys.stderr.write("### Finish extract_non_redundance_AA_seq.py ####\n\n\n")

# end: main

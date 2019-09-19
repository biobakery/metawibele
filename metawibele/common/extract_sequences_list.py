#!/usr/bin/env python
##########################################################################
# Function: Extract gene sequence based on sequnce ID 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Version: Version 1.0 09/24/2018
##########################################################################
import sys
import os
import os.path
import re
import argparse


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Extract gene sequence based on clusterID 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-l', help='input sequence ID info file', required=True)
	parser.add_argument('-i', help='input sequence file', required=True)
	parser.add_argument('-o', help='output sequence file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect sequence info 
#==============================================================
def collect_seq_info (list_file, seq_file, outfile):
	ids = {}
	open_file = open(list_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		ids[info[0]] = ""
	# foreach line
	open_file.close()

	seqs = {}
	open_file = open(seq_file, "r")
	myid = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([\S]+)", line)
			myid = mym.group(1)
			if myid in ids:
				if not myid in seqs:
					seqs[myid] = ""
			continue
		if myid in seqs:
			seqs[myid] = seqs[myid] + line
	# foreach line
	open_file.close()
	
	# output file:
	open_out = open(outfile, "w")
	for myid in sorted(ids.keys()):
		if myid in seqs:
			open_out.write(">" + myid + "\n" + seqs[myid] + "\n")
		else:
			print("No sequence!\t" + myid)
	open_out.close()
# collect_seq_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start extract_sequences.list.py -l " + values.l + " ####\n")
	collect_seq_info (values.l, values.i, values.o)
	sys.stderr.write("\n### Finish extract_sequences.list.py ####\n\n")

# end: main

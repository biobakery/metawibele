#!/usr/bin/env python

"""
MetaWIBELE: extract_non_redundance_seq module
Extract representative sequences from the clusters based on representative IDs

Copyright (c) 2019 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
import os
import re
import argparse

try:
	from metawibele import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


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
			config.logger.info ("WARNING! No sequences for gene_id: " + myid)
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-r', help='input reference sequences file', required=True)
	parser.add_argument('-i', help='input total sequences file', required=True)
	parser.add_argument('-o', help='output specific sequences file', required=True)
	values=parser.parse_args()


	config.logger.info ("### Start extract_non_redundance_seq step ####")
	
	### collect cluster id info ###
	config.logger.info ("Get clustering info ......starting")
	seqid = collect_seq_id (values.r)
	config.logger.info ("Get clustering info ......done")
	
	### Output non-redundance protein sequence
	config.logger.info ("Output sequence info ......starting")
	output_info (values.i, seqid, values.o)
	config.logger.info ("Output sequence info ......done")

	config.logger.info ("### Finish extract_non_redundance_seq step ####")

# end: main

if __name__ == '__main__':
	main()

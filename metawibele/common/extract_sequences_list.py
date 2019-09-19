#!/usr/bin/env python

"""
MetaWIBELE: extract_sequences_list module
Extract gene sequence based on sequence ID

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

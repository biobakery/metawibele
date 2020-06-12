#!/usr/bin/env python

"""
MeteWIBELE: finalize prioritization module

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

try:
	from metawibele import config
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Finalize the annotation results 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--input",
	                    help='input prioritization file',
						required = True)
	parser.add_argument('-o', "--output",
	                    help='output formated file',
	                    required=True)
	values = parser.parse_args()

	return values
# get_args


#==============================================================
# collect prioritization info
#==============================================================
def finalize_prioritization (input_file, output_file):
	title = "TID\tfamilyID\tevidence\tvalue\trank\tdescription\tnote"
	open_out = open(output_file, "w")
	open_out.write(title + "\n")

	open_file = open(input_file, "r")
	titles = {}
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[info.index(item)] = item
	mynum = 0
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		tmp = myid.split("|")
		mypid = tmp[0]
		mynote = "NA"
		if len(tmp) > 1:
			mynote = tmp[1]
		myindex = 1
		values = {}
		ranks = {}
		while myindex < len(info) -1:
			myname = titles[myindex]
			if re.search("__value$", myname):
				myname = re.sub("__value$", "", myname)
				values[myname] = info[myindex]
			if re.search("__percentile$", myname):
				myname = re.sub("__percentile$", "", myname)
				ranks[myname] = info[myindex]
			myindex = myindex + 1
		meta_name = titles[myindex]
		meta_rank = info[myindex]
		for mytype in sorted(values.keys()):
			mynum = mynum + 1
			mystr = str(mynum) + "\t" + mypid + "\t" + mytype + "\t" + values[mytype]
			myrank = "NaN"
			if mytype in ranks:
				myrank = ranks[mytype]
			mystr = mystr + "\t" + myrank + "\tranking based on single evidence"
			if mynote == "NA":
				mystr = mystr + "\t"
			else:
				mystr = mystr + "\t" + mynote
			open_out.write(mystr + "\n")
		mynum = mynum + 1
		mystr = str(mynum) + "\t" + mypid + "\t" + meta_name + "\t" + meta_rank + "\t" + meta_rank + "\t" + "meta ranking based on multiple evidences"
		if mynote == "NA":
			mystr = mystr + "\t"
		else:
			mystr = mystr + "\t" + mynote
		open_out.write(mystr + "\n")

	# foreach line
	open_file.close()
	open_out.close()
# function finalize_prioritization


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start finalize_prioritization.py -i " + values.input + " ####\n")


	### finalize prioritization for peptide families ###
	sys.stderr.write("\nFinalize prioritization for protein families ......starting\n")
	finalize_prioritization (values.input, values.output)
	sys.stderr.write("\nFinalize prioritization for protein families ......done\n")

	sys.stderr.write("### Finish finalize_prioritization.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

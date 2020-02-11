#!/usr/bin/env python

"""
MetaWIBELE: extract_abundance_feature_subset module
Extract the subset table based on specific feature

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
import math

from metawibele import utilities

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', help='input entire file', required=True)
	parser.add_argument('-s', help='input feature file', required=True)
	parser.add_argument('-o', help='output subset file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# extract subset 
#==============================================================
def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data)/n # in Python 2 use sum(data)/float(n)

def extract_subset_info (raw_file, feature_file, outfile):	
	# collect feature info
	feature = {}
	open_file = open(feature_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^" + utilities.PROTEIN_FAMILY_ID, line) or re.search("^" + utilities.PROTEIN_ID, line):
			continue
		info = line.split("\t")
		feature[info[0]] = ""
	# foreach line
	open_file.close()

	titles = {}
	mynum = 0
	sample_flt = {}
	open_file = open(raw_file, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	open_out.write(line)
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		tmp = myid.split("|")
		if tmp[0] in feature:
			open_out.write(line + "\n")
	# foreach file
	open_file.close()
	open_out.close()

# extract_subset_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start extract_abundance_feature_subset.py -i " + values.i + " ####\n")
	
	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	extract_subset_info (values.i, values.s, values.o)
	sys.stderr.write("Get abundance info ......done\n")
	

	sys.stderr.write("### Finish extract_abundance_feature_subset.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

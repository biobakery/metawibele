#!/usr/bin/env python

"""
MetaWIBELE: extract_abundance_sample_subset module
Extract the subset of abundance table based on specific samples

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

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', help='input abundance table file', required=True)
	parser.add_argument('-s', help='input sample file', required=True)
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

def extract_subset_info (raw_file, sample_file, outfile):	
	# collect sample info
	samples = {}
	open_file = open(sample_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		samples[line] = ""
	# foreach line
	open_file.close()

	titles = {}
	mynum = 0
	sample_flt = {}
	open_file = open(raw_file, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	mytitle = info[0]
	titles[0] = info[0]
	myindex = 1
	while myindex < len(info):
		item = info[myindex]
		myname = re.sub("_Abundance[\S]+", "", item)
		if myname in samples:
			titles[myindex] = item
			mytitle = mytitle + "\t" + item
			mynum = mynum + 1
			sample_flt[myname] = ""
		myindex = myindex + 1
	open_out.write(mytitle + "\n")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		mystr = info[0]
		myindex = 1
		tmp = []
		while myindex < len(info):
			if myindex in titles:
				mystr = mystr + "\t" + info[myindex]
				myvalue = info[myindex]
				if myvalue != "NA" and myvalue != "NaN" and myvalue != "nan":
					myvalue = float(myvalue)
					tmp.append(myvalue)
			myindex = myindex + 1
		# foreach sample
		if mean(tmp) == 0:
			continue
		open_out.write(mystr + "\n")	
	# foreach file
	open_file.close()
	open_out.close()

	# debug
	for mys in sorted(sample_flt.keys()):
		print(mys)
# extract_subset_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start extract_abundance_sample_subset.py -i " + values.i + " ####\n")
	
	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	extract_subset_info (values.i, values.s, values.o)
	sys.stderr.write("Get abundance info ......done\n")
	

	sys.stderr.write("### Finish extract_abundance_sample_subset.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

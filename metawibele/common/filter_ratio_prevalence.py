#!/usr/bin/env python

"""
MetaWIBELE: filter_ratio_prevalence module
Filtering features based on the ratios' prevalence 

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

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Filtering features based on prevalence 
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', 
						help='input raw ratio file', 
						required=True)
	parser.add_argument('-b', 
						help='input tranformed ratio file', 
						required=True)
	parser.add_argument('-f', 
						help='specify whether to do filtering based on prevalence: no prevalence filtering[no], filter out prevalence < 0.10[0.10]; default=[0.10]', 
						required=True,
						default="0.10")	
	parser.add_argument('-r', 
						help='whether remove infinite values', 
						choices=["no", "yes"], 
						default="no")
	parser.add_argument('-o', 
						help='output abundance file', 
						required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# filtering features
#==============================================================
def filter_feature (prevalence_flt, inf_flt, rawfile, transfile, outfile):
	open_file = open(rawfile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	sample_num = len(info) - 1
	titles = {}
	features = {}
	myindex = 0
	while myindex < len(info):
		item = info[myindex]
		titles[myindex] = item
		myindex = myindex + 1
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myindex = 1
		mynum = 0
		while myindex < len(info):
			myabu = info[myindex]
			if inf_flt == "yes":
				myabu = re.sub("Inf", "NaN", myabu)
				myabu = re.sub("inf", "NaN", myabu)
			if myabu != "NA" and myabu != "NaN" and myabu != "nan":
				myabu = float(info[myindex])
				if myabu:
					mynum = mynum + 1
			myindex = myindex + 1
		# foreach sample
		if prevalence_flt != "no":
			mypre = float(mynum) / float(sample_num)
			if mypre < float(prevalence_flt):
				continue
		features[myid] = ""
	# foreach line
	open_file.close()
	
	open_file = open(transfile, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	open_out.write(line)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if not myid in features:
			continue
		if inf_flt == "yes":
			line = re.sub("Inf", "NaN", line)
			line = re.sub("inf", "NaN", line)
		open_out.write(line + "\n")
	# foreach line
	open_file.close()
	open_out.close()

# filter_feature


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start filter_ratio_prevalence.py -a " + values.a + " ####\n")
	
	
	### filter info ###
	sys.stderr.write("Filter info ......starting\n")
	filter_feature (values.f, values.r, values.a, values.b, values.o)
	sys.stderr.write("Filter info ......done\n")


	sys.stderr.write("### Finish filter_ratio_prevalence.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

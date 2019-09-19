#!/usr/bin/env python

"""
MetaWIBELE: abundance_filtering module
Filter out potential bad proteins

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

try:
	from metawibele import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Filter out potential bad proteins
"""

def get_args (): 	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', help='input annotation file for all proteins', required=True)
	parser.add_argument('-f', help='specify whether to filter based on annotation, e.g. no | good', default="good")
	parser.add_argument('-a', help='input abundance table', required=True)
	parser.add_argument('-o', help='output abudance info for specified protein families', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect annotation info 
#==============================================================
def collect_annotation_info (ann_file, spe_note):	
	ann_cluster = {}
	titles = {}
	open_file = open(ann_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		mynote = info[titles["note"]]
		if spe_note != "no":
			if not re.search(mynote, spe_note):
				continue
		ann_cluster[myid] = ""
	# foreach line
	open_file.close()
	return ann_cluster
# collect_annotation_info


#==============================================================
# collect abundance info 
#==============================================================
def collect_abundance_info (cluster, abun_file, outfile):	
	# collect abundance info
	sys.stderr.write("Collect abundance info ......\n")
	open_file = open(abun_file, "r")
	open_out = open(outfile, "w")
	title = open_file.readline()
	open_out.write(title)
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if not info[0] in cluster:
			continue
		open_out.write(line + "\n")	
	# foreach file
	open_file.close()
	open_out.close()
# collect_abundance_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start abundance_filtering.py -a " + values.a + " ####\n")
	
	### collect stat abundance info ###
	sys.stderr.write("Get bioactivity info ......starting\n")
	cluster = collect_annotation_info (values.i, values.f)
	sys.stderr.write("Get bioactivity info ......done\n")

	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	collect_abundance_info (cluster, values.a, values.o)
	sys.stderr.write("Get abundance info ......done\n")
	

	sys.stderr.write("### Finish abundance_filtering.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python
##########################################################################
# Function: Extract the subset table based on specific feature 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 06/27/2019
##########################################################################
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
		if myid in feature:
			open_out.write(line + "\n")
	# foreach file
	open_file.close()
	open_out.close()

# extract_subset_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start extract_abundance_feature_subset.py -i " + values.i + " ####\n")
	
	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	extract_subset_info (values.i, values.s, values.o)
	sys.stderr.write("Get abundance info ......done\n")
	

	sys.stderr.write("### Finish extract_abundance_feature_subset.py ####\n\n\n")

# end: main

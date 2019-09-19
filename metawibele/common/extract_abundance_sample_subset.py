#!/usr/bin/env python
##########################################################################
# Function: Extract the subset of abundance table based on specific samples 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 06/10/2019
##########################################################################
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
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start extract_abundance_sample_subset.py -i " + values.i + " ####\n")
	
	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	extract_subset_info (values.i, values.s, values.o)
	sys.stderr.write("Get abundance info ......done\n")
	

	sys.stderr.write("### Finish extract_abundance_sample_subset.py ####\n\n\n")

# end: main

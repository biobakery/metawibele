#!/usr/bin/env python
##########################################################################
# Function: Prepare heatmap input
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Version: Version 1.0
##########################################################################
import sys
import os
import os.path
import re
import argparse
import statistics

from metawibele import config
from metawibele import utilities

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Prepare heatmap input
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', help='input normalized abundance file', required=True)
	parser.add_argument('-o', help='output formated abundance table', required=True)    
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect abundance info
#==============================================================
def collect_abundance_info (abu_file, sample_info, outfile): 
	samples = {}
	abundance = {}
	open_file = open(abu_file, "r")
	titles = {}
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[info.index(item)] = item
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myindex = 1
		while myindex < len(info):
			mys = titles[myindex]
			if mys in sample_info:
				mys = sample_info[mys] + "__" + mys
			mys = re.sub("nonIBD.dysbiosis", "nonIBD", mys)
			mys = re.sub("nonIBD.non_dysbiosis", "nonIBD", mys)
			samples[mys] = ""
			if not myid in abundance:
				abundance[myid] = {}
			abundance[myid][mys] = info[myindex]
			myindex = myindex + 1
		# foreach sample
    # foreach line
	open_file.close()
	
	# output info
	outfile1 = re.sub(".tsv", ".condition.tsv", outfile)
	open_out = open(outfile, "w")
	open_out1 = open(outfile1, "w")
	mystr = ""
	title = "familyID"
	open_out1.write("Group\tCondition\n")
	for mys in sorted(samples.keys()):
		mym = re.search("^([\S]+)__", mys)
		mycon = mym.group(1)
		#if mystr == "":
		#	mystr = mycon
		#else:
		#	mystr = mystr + "\t" + mycon
		open_out1.write(mys + "\t" + mycon + "\n")
		title = title + "\t" + mys
	open_out.write(title + "\n")
	#open_out1.write(mystr + "\n")
	open_out1.close()
	for myid in sorted(abundance.keys()):
		mystr = myid 
		for mys in sorted(samples.keys()):
			if mys in abundance[myid]:
				mystr = mystr + "\t" + str(abundance[myid][mys])
			else:
				mystr = mystr + "\t0"
		# foreach sample
		open_out.write(mystr + "\n")
	# foreach 
	open_out.close()
# function collect_cluster_info


###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start prepare_abundance_heatmap.py -a " + values.a + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	sample_info = utilities.sample_info (config.metadata, config.study)
	collect_abundance_info (values.a, sample_info, values.o)
	sys.stderr.write("Get info ......done\n")


	sys.stderr.write("### Finish prepare_abundance_heatmap.py ####\n\n\n")

# end: main

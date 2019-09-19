#!/usr/bin/env python

"""
MetaWIBELE: prepare_abundance_heatmap module
Prepare heatmap input

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

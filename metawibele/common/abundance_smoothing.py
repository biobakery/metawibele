#!/usr/bin/env python

"""
MetaWIBELE: abundance_smoothing module
Zero values were additively smoothed by half the smallest non-zero measurement on a per-sample basis and log transform

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
Zero values were additively smoothed by half the smallest non-zero measurement on a per-sample basis and log transform
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', 
						help='input abundance file', 
						required=True)
	parser.add_argument('-t', 
						help='specify the method for smoothing: the same, small smoothing factor across samples [fixed], half the smallest non-zero value per sample[unfixed]; default=[fixed]', 
						choices=["fixed", "unfixed"],
						default="fixed")
	parser.add_argument('-f', 
						help='specify whether to do filtering based on prevalence: no prevalence filtering[no], filter out prevalence < 0.10[0.10]; default=[0.10]', 
						required=True, 
						default="0.10")
	parser.add_argument('-o', 
						help='output smoothed abundance file', 
						required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster abundance
#==============================================================
def collect_cluster_abundance (abufile):
	samples = {}
	samples_tmp = {}
	titles = {}
	open_file = open(abufile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 1
	while myindex < len(info):
		titles[myindex] = info[myindex]
		myindex = myindex + 1
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myindex = 1
		while myindex < len(info):
			mys = titles[myindex]
			if float(info[myindex]) == 0:
				myindex = myindex + 1
				continue
			if not mys in samples_tmp:
				samples_tmp[mys] = []
			samples_tmp[mys].append(float(info[myindex]))
			myindex = myindex + 1
	# foreach line
	open_file.close()

	# the smallest non-zero abundance per sample
	mins = []
	for mys in samples_tmp:
		mymin = min(samples_tmp[mys])
		samples[mys] = mymin * 0.5
		mins.append(mymin)
	# foreach sample

	return samples, min(mins)
#collect_cluster_abundance


#==============================================================
# smooth abundance
#==============================================================
def smooth_abundance (smooth_method, prevalence_flt, abufile, samples, min_smooth, outfile):
	values = {}
	diagnosis = {}
	types = {}
	outfile1 = re.sub(".tsv", ".refined.tsv", abufile)
	outfile2 = re.sub(".tsv", ".log.tsv", outfile)
	open_file = open(abufile, "r")
	open_out = open(outfile, "w")
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile2, "w")
	line = open_file.readline()
	open_out.write(line)
	open_out1.write(line)
	open_out2.write(line)
	line = line.strip()
	info = line.split("\t")
	sample_num = len(info) - 1
	titles = {}
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
		mystr = info[0]
		mystr2 = info[0]
		myindex = 1
		mynum = 0
		while myindex < len(info):
			myabu = info[myindex]
			myabu2 = info[myindex]
			if myabu != "NA" and myabu != "NaN" and myabu != "nan":
				myabu = float(info[myindex])
				mys = titles[myindex]
				if myabu > float(config.abundance_detection_level):
					mynum = mynum + 1
				else:
					if smooth_method == "fixed":
						myabu = float(min_smooth)
					else:
						if mys in samples:
							myabu = float(samples[mys])
						else:
							# debug
							print("No sample smoothing value!\t" + info[0] + "\t" + mys)
							myabu = float("NaN")
				if not math.isnan(myabu):
					myabu2 = math.log(myabu)
			mystr = mystr + "\t" + str(myabu)
			mystr2 = mystr2 + "\t" + str(myabu2)
			myindex = myindex + 1
		# foreach sample
		if prevalence_flt != "no":
			mypre = float(mynum) / float(sample_num)
			if mypre < float(prevalence_flt):
				continue
		open_out.write(mystr + "\n")
		open_out2.write(mystr2 + "\n")
		open_out1.write(line + "\n")
	# foreach line
	open_file.close()
	open_out.close()
	open_out1.close()
	open_out2.close()

# smooth_abundance


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start abundance_smoothing.py -i " + values.i + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	samples, min_smooth = collect_cluster_abundance (values.i)
	sys.stderr.write("Get info ......done\n")
	
	### convert info ###
	sys.stderr.write("Convert info ......starting\n")
	smooth_abundance (values.t, values.f, values.i, samples, min_smooth, values.o)
	sys.stderr.write("Output info ......done\n")


	sys.stderr.write("### Finish abundance_smoothing.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

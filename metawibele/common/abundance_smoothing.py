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
						choices=["fixed", "sample", "feature"],
						default="feature")
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
	'''
	:param abufile: abundance file
	:return: abundance table
	'''

	samples = {}
	features = {}
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
		myid = info[0]
		myindex = 1
		feature_tmp = []
		while myindex < len(info):
			mys = titles[myindex]
			myabu = info[myindex]
			if myabu != "NA" and myabu != "NaN" and myabu != "nan":
				myabu = float(info[myindex])
				if myabu > float(config.abundance_detection_level):
					if not mys in samples_tmp:
						samples_tmp[mys] = []
					samples_tmp[mys].append(myabu)
					feature_tmp.append(myabu)
			myindex = myindex + 1
		# foreach sample
		if len(feature_tmp) > 0:
			features[myid] = min(feature_tmp) * 0.5
	# foreach line
	open_file.close()

	# the smallest non-zero abundance per sample
	mins = []
	for mys in samples_tmp:
		mymin = min(samples_tmp[mys])
		samples[mys] = mymin * 0.5
		mins.append(mymin)
	# foreach sample
	fixed_min = min(mins)

	return  fixed_min, samples, features
#collect_cluster_abundance


#==============================================================
# smooth abundance
#==============================================================
def smooth_abundance (smooth_method, prevalence_flt, abufile, fixed_min, samples, features, outfile):
	'''
	:param smooth_method: smoothing method [fixed, sample, feature]
	:param prevalence_flt: thredshold for prevalence filtering
	:param abufile: abundance file
	:param fixed_min: global fixed pseudocount
	:param samples: per-sample pseudocounts
	:param features: per-feature pseudocount
	:param outfile: smoothed file name
	:return: smoothed file
	'''
	
	myname = os.path.basename(abufile)
	mydir = os.path.dirname(outfile)
	outfile1 = re.sub(".tsv", ".refined.tsv", myname)
	outfile1 = os.path.join(mydir, outfile1)
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
		myid = info[0]
		mystr = myid
		mystr2 = myid
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
						myabu = float(fixed_min)
					if smooth_method == "sample":
						if mys in samples:
							myabu = float(samples[mys])
						else:
							# debug
							print("Per-sample pseudocount is unavailable!\t" + myid + "\t" + mys)
							myabu = float("NaN")
					if smooth_method == "feature":
						if myid in features:
							myabu = float(features[myid])
						else:
							# debug
							print("Per-feature pseudocount is unavailable!\t" + myid)
							myabu = float("NaN")
				if not math.isnan(myabu):
					myabu2 = math.log(myabu)
			mystr = mystr + "\t" + str(myabu)
			mystr2 = mystr2 + "\t" + str(myabu2)
			myindex = myindex + 1
		# foreach sample

		# if prevalence_flt != "no"
		try:
			prevalence_flt = float(prevalence_flt)
			mypre = float(mynum) / float(sample_num)
			if mypre < prevalence_flt:
				continue
		except:
			mypre = float(mynum) / float(sample_num)
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
	fixed_min, samples, features = collect_cluster_abundance (values.i)
	sys.stderr.write("Get info ......done\n")
	
	### smooth info ###
	sys.stderr.write("Smooth info ......starting\n")
	smooth_abundance (values.t, values.f, values.i, fixed_min, samples, features, values.o)
	sys.stderr.write("Smooth info ......done\n")


	sys.stderr.write("### Finish abundance_smoothing.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

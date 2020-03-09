#!/usr/bin/env python3

"""
MetaWIBELE: rna_dna_ratio module
Get the ratio value of RNA and DNA 

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
import statistics
import numpy as np

try:
	from metawibele import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
			" Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Get the ratio value of RNA and DNA 
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-d', 
						help='input relative abundance file for DNA', 
						required=True)
	parser.add_argument('-r', 
						help='input relative abundance file for RNA', 
						required=True)
	parser.add_argument('-o', 
						help='output ratio table', 
						required=True)    
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# basis stattistics info
#==============================================================
def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data)/n # in Python 2 use sum(data)/float(n)

def _ss(data):
	"""Return sum of square deviations of sequence data."""
	c = mean(data)
	ss = sum((x-c)**2 for x in data)
	return ss

def stddev(data, ddof=0):
	"""Calculates the population standard deviation
	by default; specify ddof=1 to compute the sample
	standard deviation."""
	n = len(data)
	if n < 2:
		raise ValueError('variance requires at least two data points')
	ss = _ss(data)
	pvar = ss/(n-ddof)
	return pvar**0.5


#==============================================================
# collect cluster abundance
#==============================================================
def collect_cluster_abundance (abufile):
	'''
	:param abufile: abundance file
	:return: abundance table
	'''

	abundance = {}
	titles = {}
	ids = []
	smooth = {}
	smooth_tmp = {}
	open_file = open(abufile, "r")
	title = open_file.readline()
	title = title.strip()
	info = title.split("\t")
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
		ids.append(myid)
		myindex = 1
		while myindex < len(info):
			mys = titles[myindex]
			if not myid in abundance:
				abundance[myid] = {}
			abundance[myid][mys] = info[myindex]
			if info[myindex] != "NaN" and info[myindex] != "NA" and info[myindex] != "nan":
				if float(info[myindex]) > float(config.abundance_detection_level):
					if not mys in smooth_tmp:
						smooth_tmp[mys] = []
					smooth_tmp[mys].append(float(info[myindex]))
			myindex = myindex + 1
	# foreach line
	open_file.close()

	# find the smallest value for each sample
	mymin = []
	mymax = []
	for mys in smooth_tmp.keys():
		myvalue = min(smooth_tmp[mys]) / 2
		smooth[mys] = myvalue
		mymin.append(myvalue)
		myvalue = max(smooth_tmp[mys])
		mymax.append(myvalue)
	# foreach sample
	return title, ids, abundance, smooth, min(mymin), max(mymax)
#collect_cluster_abundance


#==============================================================
# calculate ratio value and smooth with pseudocount
#==============================================================
def calculate_ratio (abundance_dna, rna_file, outfile):
	'''
	:param abundance_dna: DNA abundance desc
	:param rna_file: RNA abundance file
	:outfile: ratio file name
	:return: RNA-DNA ratio file and corresponding log transformed files
	'''

	samples = {}
	infs = {}
	zeros = {}
	open_file = open(rna_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 1
	while myindex < len(info):
		mys = info[myindex]
		samples[myindex] = mys
		myindex = myindex + 1
	# foreach sample
	outfile2 = re.sub(".tsv", ".log.tsv", outfile)
	open_out = open(outfile, "w")
	open_out2 = open(outfile2, "w")
	open_out.write(line + "\n")
	open_out2.write(line + "\n")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if not myid in abundance_dna:
			continue
		mystr = myid
		mystr2 = myid
		myindex = 1
		myna = 0
		while myindex < len(info):
			myrna = info[myindex]
			mys = samples[myindex]
			if not mys in abundance_dna[myid]:
				# debug
				print("No DNA info!\t" + myid + "\t" + mys)
				myvalue = "NaN"
				mystr = mystr + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
				myna = myna + 1
				myindex = myindex + 1
				continue
			mydna = abundance_dna[myid][mys]
			if myrna == "NA" or myrna == "NaN" or myrna == "nan":
				myvalue = "NaN"
				mystr = mystr + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
				myna = myna + 1
				myindex = myindex + 1
				continue
			if mydna == "NA" or mydna == "NaN" or mydna == "nan":
				myvalue = "NaN"
				mystr = mystr + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
				myna = myna + 1
				myindex = myindex + 1
				continue
			myrna = float(myrna)
			mydna = float(mydna)
			if myrna == 0 and mydna == 0:
				myvalue = "NaN"
				mystr = mystr + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
				myna = myna + 1
			if myrna == 0 and mydna != 0:
				mystr = mystr + "\t0" 
				mystr2 = mystr2 + "\tNaN"
				myna = myna + 1
				zeros[mys + "\t" + myid] = mydna
			if myrna != 0 and mydna == 0:
				myvalue = "inf"
				mystr = mystr + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
				infs[mys + "\t" + myid] = myrna
			if myrna != 0 and mydna != 0:
				myvalue = myrna / mydna
				mystr = mystr + "\t" + str(myvalue)
				mystr2 = mystr2 + "\t" +  str(math.log(myvalue))
			myindex = myindex + 1
		# foreach sample
		if myna == len(samples.keys()):
			continue
		open_out.write(mystr + "\n")
		open_out2.write(mystr2 + "\n")
	# foreach line
	open_file.close()
	open_out.close()
	open_out2.close()

	return infs, zeros
# calculate_ratio


def calculate_smoothing_value (abu_file):
	'''
	:param abu_file:
	:return: smoothed values
	'''

	# collect info
	samples = {}
	samples_tmp = {}
	features_min = {}
	features_max = {}
	features_quantile = {}
	titles = {}
	open_file = open(abu_file, "r")
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
			if myabu != "NA" and myabu != "NaN" and myabu != "nan" and myabu != "inf":
				myabu = float(info[myindex])
				if myabu != 0:
					if not mys in samples_tmp:
						samples_tmp[mys] = []
					samples_tmp[mys].append(myabu)
					feature_tmp.append(myabu)
			myindex = myindex + 1
		# foreach sample
		if len(feature_tmp) > 0:
			features_min[myid] = min(feature_tmp) * 0.5
			features_max[myid] = max(feature_tmp) * 0.99
			feature_tmp.sort()
			#features_quantile[myid] = statistics.median(feature_tmp)
			features_quantile[myid] = np.percentile(feature_tmp, 10)
	# foreach line
	open_file.close()

	# the smallest non-zero ratio per sample
	mins = []
	for mys in samples_tmp:
		mymin = min(samples_tmp[mys])
		samples[mys] = mymin * 0.5
		mins.append(mymin)
	# foreach sample
	fixed_min = min(mins)

	return features_min, features_max, features_quantile

# calculate_smoothing_value


def smooth_zero (dna_features_min, dna_features_quantile, rna_features_min, rna_features_quantile, infs, zeros, ratio_file):
	outfile1 = re.sub(".tsv", ".smooth.tsv", ratio_file)
	outfile2 = re.sub(".tsv", ".smooth.log.tsv", ratio_file)
	open_file = open(ratio_file, "r")
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile2, "w")
	line = open_file.readline()
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
		mystr1 = myid
		mystr2 = myid
		myindex = 1
		mynum = 0
		while myindex < len(info):
			myabu = info[myindex]
			myabu2 = info[myindex]
			if myabu != "NA" and myabu != "NaN" and myabu != "nan":
				myabu = float(info[myindex])
				mys = titles[myindex]
				if myabu > 0:
					mynum = mynum + 1
					if str(myabu) == "inf":
						myrna = "NA"
						mydna = "NA"
						mytmp = mys + "\t" + myid
						if mytmp in infs:
							myrna = infs[mytmp]
							#if myid in rna_features_quantile: # skip extremely small values for both DNA and RNA
							#	if float(myrna) < float(rna_features_quantile[myid]):
							#		myrna = "NA"
							#else:
							#	myrna = "NA"
						if myid in dna_features_min:
							mydna = dna_features_min[myid]
						if myrna != "NA" and mydna != "NA":
							myabu = float(myrna) / float(mydna)
						else:
							myabu = "NaN"
							myabu2 = "NaN"
				else:
					# rna == 0
					myrna = "NA"
					mydna = "NA"
					mytmp = mys + "\t" + myid
					if mytmp in zeros:
						mydna = zeros[mytmp]
						#if myid in dna_features_quantile: # skip extremely small values for both DNA and RNA
						#	if float(mydna) < float(dna_features_quantile[myid]):
						#		mydna = "NA"
						#else:
						#	mydna = "NA"
					if myid in rna_features_min:
						myrna = rna_features_min[myid]
					if myrna != "NA" and mydna != "NA":
						myabu = float(myrna) / float(mydna)
					else:
						myabu = "NaN"
						myabu2 = "NaN"
				if not math.isnan(float(myabu)):
					if float(myabu) != 0:
						myabu2 = math.log(myabu)
			mystr1 = mystr1 + "\t" + str(myabu)
			mystr2 = mystr2 + "\t" + str(myabu2)
			myindex = myindex + 1
		# foreach sample
		open_out1.write(mystr1 + "\n")
		open_out2.write(mystr2 + "\n")
	# foreach line
	open_file.close()
	open_out1.close()
	open_out2.close()

# smooth_zero


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start rna_dna_ratio_percentile.py -d " + values.d + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	title_dna, ids_dna, abundance_dna, smooth_dna, min_dna, max_dna = collect_cluster_abundance (values.d)
	infs, zeros = calculate_ratio (abundance_dna, values.r, values.o)
	dna_features_min, dna_features_max, dna_features_quantile = calculate_smoothing_value (values.d)
	rna_features_min, rna_features_max, rna_features_quantile = calculate_smoothing_value (values.r)
	#ratio_features_min, ratio_features_max, ratio_features_quantile = calculate_smoothing_value (values.o)
	sys.stderr.write("Get info ......done\n")
	
	### smooth zero ###
	sys.stderr.write("Smooth info ......starting\n")
	smooth_zero (dna_features_min, dna_features_quantile, rna_features_min, rna_features_quantile, infs, zeros, values.o)
	sys.stderr.write("Smooth info ......done\n")


	sys.stderr.write("### Finish rna_dna_ratio_percentile.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

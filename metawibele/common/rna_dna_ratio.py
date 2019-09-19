#!/usr/bin/env python

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
	parser.add_argument('-t', 
						help='specify the method for smoothing: the same, small smoothing factor across samples [fixed], half the smallest non-zero value per sample[unfixed]; default=[fixed]', 
						choices=["fixed", "unfixed"],
						default="fixed")	
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
# calculate ratio value
#==============================================================
def convert_abundance_info (abundance_dna, smooth_dna, smooth_rna, min_rna, rna_file, smooth_flag, outfile): 
	samples = {}
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
	outfile1 = re.sub(".tsv", ".smooth.log.tsv", outfile)
	outfile2 = re.sub(".tsv", ".log.tsv", outfile)
	open_out = open(outfile, "w")
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile2, "w")
	open_out.write(line + "\n")
	open_out1.write(line + "\n")
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
		mystr1 = myid
		mystr2 = myid
		myindex = 1
		myna = 0
		while myindex < len(info):
			myrna = float(info[myindex])
			mys = samples[myindex]
			if not mys in abundance_dna[myid]:
				# debug
				print("No DNA info!\t" + myid + "\t" + mys)
				myvalue = "NaN"
				mystr = mystr + "\t" + myvalue
				mystr1 = mystr1 + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
				myna = myna + 1
				continue
			mydna = float(abundance_dna[myid][mys])
			if myrna == 0 and mydna == 0:
				myvalue = "NaN"
				mystr = mystr + "\t" + myvalue
				mystr1 = mystr1 + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
				myna = myna + 1
			if myrna == 0 and mydna != 0:
				mystr = mystr + "\t0" 
				mystr2 = mystr2 + "\tNaN"
				if smooth_flag == "fixed":
					myvalue = math.log(float(min_rna) / mydna)
					#myvalue = math.log(float(min_rna))
				else:
					if mys in smooth_rna:
						myvalue = math.log(smooth_rna[mys] / mydna)
					else:
						# debug
						print("No smoothed RNA info!\t" + mys)
						myvalue = "NaN"
				mystr1 = mystr1 + "\t" + str(myvalue)
				myna = myna + 1
			if myrna != 0 and mydna == 0:
				#if mys in smooth_dna:
				#	mysmall = smooth_dna[mys]
				#	myvalue = myrna / mysmall
				#	if trans_flag != "no":
				#		if trans_flag == "log":
				#			myvalue = math.log(myvalue)
				#else:
				myvalue = "inf"
				mystr = mystr + "\t" + myvalue
				mystr1 = mystr1 + "\t" + myvalue
				mystr2 = mystr2 + "\t" + myvalue
			if myrna != 0 and mydna != 0:
				myvalue = myrna / mydna
				mystr = mystr + "\t" + str(myvalue)
				mystr1 = mystr1 + "\t" +  str(math.log(myvalue))
				mystr2 = mystr2 + "\t" +  str(math.log(myvalue))
			myindex = myindex + 1
		# foreach sample
		if myna == len(samples.keys()):
			continue
		open_out.write(mystr + "\n")
		open_out1.write(mystr1 + "\n")
		open_out2.write(mystr2 + "\n")
	# foreach line
	open_out.close()
	open_out1.close()
	open_out2.close()
	open_file.close()
# convert_abundance_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start rna_dna_ratio.py -d " + values.d + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	title_dna, ids_dna, abundance_dna, smooth_dna, min_dna, max_dna = collect_cluster_abundance (values.d)
	title_rna, ids_rna, abundance_rna, smooth_rna, min_rna, max_rna = collect_cluster_abundance (values.r)
	abundance_rna = {}
	sys.stderr.write("Get info ......done\n")
	
	### convert info ###
	sys.stderr.write("Convert info ......starting\n")
	#mysmooth = float(min_rna) / float(max_dna)
	convert_abundance_info (abundance_dna, smooth_dna, smooth_rna, min_rna, values.r, values.t, values.o)
	sys.stderr.write("Output info ......done\n")


	sys.stderr.write("### Finish rna_dna_ratio.py ####\n\n\n")

# end: main

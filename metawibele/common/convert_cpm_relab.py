#!/usr/bin/env python3

"""
MetaWIBELE: convert_cpm_relab module
Convert CPM to relative abundance [0,1]

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
						help='specify data type, e.g. log | cpm', 
						default="cpm")
	parser.add_argument('-o', 
						help='output smoothed abundance file', 
						required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# convert abundance
#==============================================================
def convert_abundance (abufile, data_type, outfile):
	'''
	:param abufile: abundance file
	:return: abundance table
	'''

	open_file = open(abufile, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	open_out.write(line)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		mystr = myid 
		myindex = 1
		while myindex < len(info):
			myabu = info[myindex]
			if myabu != "NA" and myabu != "NaN" and myabu != "nan":
				if data_type == "cpm":
					myabu = float(info[myindex]) / 1000000
				if data_type == "log":
					a = math.exp(float(myabu))
					a = a / 1000000
					myabu = math.log(a)
			mystr = mystr + "\t" + str(myabu)
			myindex = myindex + 1
		open_out.write (mystr + "\n")
	# foreach line
	open_file.close()
	open_out.close()

# convert_abundance


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start convert_cpm_relab.py -i " + values.i + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	convert_abundance (values.i, values.t, values.o)
	sys.stderr.write("Get info ......done\n")
	

	sys.stderr.write("### Finish convert_cpm_relab.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

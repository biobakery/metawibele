#!/usr/bin/env python

"""
MetaWIBELE: metadata_format module
Format the metadata info 

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

try:
	from metawibele import config
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Format the metadata info and specify reference status
"""
	
def get_args():	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', help='input meta data file', required=True)
	parser.add_argument('-o', help='output formatted meta file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# format metadata
#==============================================================
def format_metadata (infile, outfile):
	titles = {}
	if not os.path.isfile(infile):
		print ("File not exist!\t" + infile)
		return
	open_file = open(infile, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	line = re.sub("External_ID", "ID", line)
	line = re.sub("SID", "ID", line)
	open_out.write(line)
	line = line.strip()
	info = line.split("\t")
	for item in info:
		myindex = info.index(item)
		titles[item] = myindex
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		for mymeta in config.ref_status.keys():
			if mymeta in titles:
				myindex = titles[mymeta]
				myindex = titles[mymeta]
				item = info[myindex]
				tmp2 = config.ref_status[mymeta].split(",")
				myrefs = []
				for mytmp in tmp2:
					if re.search("_vs_", mytmp):
						mym = re.search("_vs_([\S]+)", mytmp)
						myrefs.append(mym.group(1))
					else:
						myrefs.append(mytmp)
				if item in myrefs:
					info[myindex] = "a_" + info[myindex]
		# foreach metadata
		mystr = "\t".join(info)
		open_out.write(mystr + "\n")
	# foreach line
	open_file.close()
	open_out.close()
# format_metadata


#==============================================================
###########  Main processing ############
#==============================================================
def main():

	### get arguments ###

	values = get_args()


	sys.stderr.write("### Start metadata_format.py -o " + values.o + " ####\n")
	
	### collect sample info ###
	sys.stderr.write("Format info ......starting\n")
	format_metadata (values.i, values.o)
	sys.stderr.write("Format info ......done\n")
	
	sys.stderr.write("### Finish metadata_format.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

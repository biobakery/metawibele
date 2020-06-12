#!/usr/bin/env python

"""
MetaWIBELE: filter_clusters module
Filter out specific clusters

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


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Filter out specific clusters 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-l', help='input cluster ID info file', required=True)
	parser.add_argument('-i', help='input raw file', required=True)
	parser.add_argument('-o', help='output refined file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect info 
#==============================================================
def collect_info (list_file, info_file, outfile):
	ids = {}
	open_file = open(list_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		ids[info[0]] = ""
	# foreach line
	open_file.close()

	open_file = open(info_file, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	open_out.write(line)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		tmp = myid.split("|")
		if tmp[0] in ids:
			continue
		open_out.write(line + "\n")
	# foreach line
	open_file.close()
	open_out.close()

# collect_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start filter_clusters.py -l " + values.l + " ####\n")
	collect_info (values.l, values.i, values.o)
	sys.stderr.write("\n### Finish filter_clusters.py ####\n\n")

# end: main

if __name__ == '__main__':
	main()

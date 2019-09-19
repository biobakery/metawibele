#!/usr/bin/env python

"""
MetaWIBELE: format_mtx module
Format MGX table for raw sequences

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
Format MGX table for raw sequences
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-s', help='input sample list', required=True)
	parser.add_argument('-b', help='input basic sequence info', required=True)
	parser.add_argument('-o', help='output formated table', required=True)    
	values = parser.parse_args()
	return values
# get_args



#==============================================================
# format squences info for MGX samples
#==============================================================
def format_mgx_info (sample_file, base_file,outfile):	
	samples = {}
	open_file = open(sample_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		samples[info[0]] = ""
	# foreach line
	open_file.close()

	open_out = open(outfile, "w")
	open_file = open(base_file, "r")
	line = open_file.readline()
	line = line.strip()
	open_out.write("SID\t" + line + "\tfile_name\n")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		for mys in sorted(samples.keys()):
			myfile = mys + ".tar"
			if re.search("_P$", mys):
				myfile = mys + ".fastq.bz2"
			open_out.write(mys + "\t" + line + "\t" + myfile + "\n")
	# foreach line
	open_out.close()
	open_file.close()

# format_mgx_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start format_mgx_table.py -s " + values.s + " ####\n")
	

	### format info ###
	sys.stderr.write("Format info ......starting\n")
	format_mgx_info (values.s, values.b, values.o)
	sys.stderr.write("Format info ......done\n")
	

	sys.stderr.write("### Finish format_mgx_table.py ####\n\n\n")

# end: main

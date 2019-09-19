#!/usr/bin/env python
##########################################################################
# Function: Format MGX table for raw sequences
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 07/07/2019
##########################################################################
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
				myfile = mys + ".fastq.gz"
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


	sys.stderr.write("### Start format_mgx.py -s " + values.s + " ####\n")
	

	### format info ###
	sys.stderr.write("Format info ......starting\n")
	format_mgx_info (values.s, values.b, values.o)
	sys.stderr.write("Format info ......done\n")
	

	sys.stderr.write("### Finish format_mgx.py ####\n\n\n")

# end: main

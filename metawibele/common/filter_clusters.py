#!/usr/bin/env python
##########################################################################
# Function: Filter out specific clusters 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 07/17/2019
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
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start filter_clusters.py -l " + values.l + " ####\n")
	collect_info (values.l, values.i, values.o)
	sys.stderr.write("\n### Finish filter_clusters.py ####\n\n")

# end: main

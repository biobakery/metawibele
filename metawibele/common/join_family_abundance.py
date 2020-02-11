#!/usr/bin/env python

"""
MetaWIBELE: join_family_abundance module
Join the abundance of families

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


#==============================================================
# join abundance info
#==============================================================
def join_abundance_info (abundance_list, outfile):	
	open_list = open(abundance_list, "r")
	samples = {}
	abundance = {}
	title = "ID"
	for myfile in open_list:
		myfile = myfile.strip()
		if not len(myfile):
			continue
		open_file = open(myfile, "r")
		titles = {}
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			if re.search("^ID", line): # title
				title = title + "\t" + "\t".join(info[1:len(info)])
				continue
			myid = info[0]
			if not myid in abundance:
				abundance[myid] = myid
			abundance[myid] = abundance[myid] + "\t" + "\t".join(info[1:len(info)])
		open_file.close()
	# foreach file	
	open_list.close()

	# output info
	open_out = open(outfile, "w")
	open_out.write(title + "\n")
	for myid in sorted(abundance.keys()):
		open_out.write(abundance[myid] + "\n")
	# foreach cluster
	open_out.close()
# function join_abundance_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-l', help='input the list of abundance info for families', required=True)
	parser.add_argument('-o', help='output joint abundance file', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start join_family_abundance.py -l " + values.l + " ####\n")
	
	
	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	join_abundance_info(values.l, values.o)
	sys.stderr.write("Get abundance info ......done\n")


	sys.stderr.write("### Finish join_family_abundance.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

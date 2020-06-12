#!/usr/bin/env python

"""
MetaWIBELE: extract_complete_ORF_seq module
Extract sequences of complete ORFs

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
# collect gene info
#==============================================================
def collect_gene_info (infile, spe_type):
	genes = {}
	titles = {}
	open_file = open(infile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		gene = info[titles["GID"]]
		mypartial = info[titles["partial"]]
		feature = info[titles["feature"]]
		if feature != "CDS":
			continue
		if mypartial == "00":	# complete genes
			flag = "complete"
		if mypartial == "10": # lose stop codon
			flag = "no_stop_codon"
		if mypartial == "01": # lose start codon
			flag = "no_start_codon"
		if mypartial == "11": # no stop and start
			flag = "no_start_stop"
		if flag == spe_type:
			genes[gene] = flag
	# foreach line
	open_file.close()
	return genes
# collect_gene_info


#==============================================================
# output info
#==============================================================
def output_info (genes, infile, outfile): 
	open_file = open(infile, "r")
	open_out = open(outfile, "w")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([\S]+)", line)
			myid = mym.group(1)
			if myid in genes:
				flag = 1
			else:
				flag = 0
		# gene id
		if flag == 1:
			open_out.write(line + "\n")
	open_file.close()
	open_out.close()
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-t', help='specify gene type, e.g. complete', required=True)
	parser.add_argument('-m', help='input sample info file', required=True)
	parser.add_argument('-i', help='input raw sequence file', required=True)
	parser.add_argument('-o', help='output subset sequence file', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start extract_complete_ORF_seq.py -i " + values.i + " ####\n")
	

	### collect gene info  ###
	sys.stderr.write("Get gene info ......starting\n")
	genes = collect_gene_info (values.m, values.t)
	sys.stderr.write("Get gene info ......done\n")
	
	### Output sequence info
	sys.stderr.write("\nOutput sequence info ......starting\n")
	output_info (genes, values.i, values.o)
	sys.stderr.write("Output sequence info ......done\n")

	sys.stderr.write("### Finish extract_complete_ORF_seq.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

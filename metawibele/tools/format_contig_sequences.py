#!/usr/bin/env python

"""
MetaWIBELE: format_contig_sequences module
Format contigs sequence info

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
import re
import argparse

from metawibele import utilities


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', help='input the path of contig files', required=True)
	parser.add_argument('-e', help='specify the file extension of contig file', default="contigs.fa")
	parser.add_argument('-o', help='output file name', required=True)
	values = parser.parse_args()
	
	return values

#==============================================================
# format contig sequences
#==============================================================
def format_contig_info (contig_path, extension, outfile):
	filelist = utilities.find_files(contig_path, extension, None)
	open_out = open(outfile, "w")
	for myfile in filelist:
		myfile = myfile.strip()
		if not len(myfile):
			continue                      
		sample = myfile
		mym = re.search("([^\/]+)$", sample)
		sample = mym.group(1)
		sample = re.sub("." + extension, "", sample)
		# collect seq info
		if not os.path.isfile(myfile):
			print("Contig file doesn't exist!\t" + myfile)
			continue
		open_contig = open(myfile, "r")
		contigs = {}
		contig_order = []
		myid = ""
		for line in open_contig:
			line = line.strip()
			if not len(line):
				continue
			if re.search("^>", line):
				mym = re.search(">([\S]+)", line)
				myid_old = mym.group(1)
				myid_new = sample + "_contig_" + mym.group(1) + "|" + sample + "|"
				myid = re.sub(myid_old + "[\s]+", myid_new, line)
				if not myid in contigs:
					contig_order.append(myid)
					contigs[myid] = ""
				continue
			contigs[myid] = contigs[myid] + line
		# foreach line
		open_contig.close()
		# output contig sequence
		for myid in contig_order:
			if myid in contigs:
				open_out.write(myid + "\n" + contigs[myid] + "\n")
		# foreach contig
	# foreach sample
	
	open_out.close()
# format_contig_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start format_contig_table.py -p " + values.p + " ####\n")
	
	### collect sequence info ###
	sys.stderr.write("Get contig and gff info ......starting\n")
	format_contig_info (values.p, values.e, values.o)
	sys.stderr.write("Get contig and gff info ......done\n")

	sys.stderr.write("### Finish format_contig_table.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

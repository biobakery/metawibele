#!/usr/bin/env python

"""
MetaWIBELE: split_family_abundance module
Split the abundance of families

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
# split abundance info
#==============================================================
def split_abundance_info (abundance_file, sample_num):	# summary_peptide_family_abundance.RPK.all.tsv
	open_file = open(abundance_file, "r")
	titles = []
	samples = {}
	abundance = {}
	myflag = "NA"
	mypos = 0
	split = 0
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^ID", line):
			myflag = "ID"
			for item in info:
				#titles[info.index(item)] = item
				titles.append(item)
			continue
		if re.search("^# Gene Family", line): # title
			myflag = "# Gene Family"
			for item in info:
				titles.append(item)
			continue
		myid = info[0]
		if not myid in abundance:
			abundance[myid] = {}
		myindex = 1
		split = 0
		while myindex < len(info):
			#if mypos % int(sample_num) == 0:	# new split
			split = split + 1
			mys = myindex
			mye = myindex + int(sample_num)
			if mye > len(info):
				mye = len(info)
			myinfo = "\t".join(info[mys:mye])
			mysample = "\t".join(titles[mys:mye])
			if not split in samples:
				samples[split] = mysample
			if not split in abundance:
				abundance[split] = {}
			abundance[split][myid] = myinfo
			myindex = mye
		# foreach sample	
	open_file.close()
	
	# output info
	for mysplit in sorted(samples.keys(), key=int):
		outfile = re.sub(".tsv", ".split" + str(mysplit) + ".tsv", abundance_file)
		open_out = open(outfile, "w")
		title = myflag + "\t" + samples[mysplit]
		open_out.write(title + "\n")
		if mysplit in abundance:
			for myid in sorted(abundance[mysplit].keys()):
				mystr = myid + "\t" + abundance[mysplit][myid]
				open_out.write(mystr + "\n")
		# foreach cluster
		open_out.close()
	# foreach split
# function split_abundance_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-a', help='input the abundance info for families', required=True)
	parser.add_argument('-n', help='the number of samples for each split file', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start split_family_abundance.py -a " + values.a + " ####\n")
	
	
	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	split_abundance_info(values.a, values.n)
	sys.stderr.write("Get abundance info ......done\n")


	sys.stderr.write("### Finish split_family_abundance.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python
##########################################################################
# Function: Split the abundance of families
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Version: Version 1.0	12/14/2018
##########################################################################
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
	mypos = 0
	split = 0
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^ID", line): # title
			for item in info:
				#titles[info.index(item)] = item
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
		title = "ID\t" + samples[mysplit]
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
if __name__ == '__main__':
	
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

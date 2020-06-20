#!/usr/bin/env python

"""
MetaWIBELE: abundance_RPK module
Normalize family abundance to RPK (read counts divided by the length of genes in kilobase to normlized for gene/protein length)

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
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Normalize family abundance to RPK (read counts divided by the length of genes in kilobase to normlized for gene/protein length)
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', help='input protein abundance (read count)', required=True)
	parser.add_argument('-l', help='input the file including the protein sequence length', required=True)
	parser.add_argument('-t', help='specify the type of length file', choices = ['fasta', 'cluster'], default='cluster')
	parser.add_argument('-o', help='output normalized abundance file', required=True)    
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_cluster_info (clust_file, file_type):	
	cluster = {}
	seqs = {}
	open_file = open(clust_file, "r")
	myclust = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			if file_type == "cluster":
				mym = re.search("cluster=([\d]+)", line)
				myclust = "Cluster_" + mym.group(1)
				mym = re.search("length=([\d]+)", line)
				mylen = mym.group(1)
				mym = re.search(">([^;]+)", line)
				myid = mym.group(1)
				if not myclust in cluster:
					cluster[myclust] = mylen
				continue
			if file_type == "fasta":
				mym = re.search(">([\S]+)", line)
				myid = mym.group(1)
				seqs[myid] = ""
		if file_type == "fasta":
			seqs[myid] = seqs[myid] + line
	# foreach line
	open_file.close()

	if file_type == "fasta":
		for myid in seqs.keys():
			mylen = len(seqs[myid])
			cluster[myid] = mylen

	return cluster
# function collect_cluster_info


#==============================================================
# normalize abundance
#==============================================================
def normalization (cluster, abun_file, outfile):	# summary_peptide_family_abundance.tsv
	open_out = open (outfile, "w")
	open_file = open (abun_file, "r")
	line = open_file.readline()
	open_out.write(line)	# output title
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myclust = info[0]
		mylen = 1
		if myclust in cluster:
			mylen = float(cluster[myclust]) * 3
			mylen = float(mylen) / 1000
			myindex = 1
			mystr = myclust
			while myindex < len(info):
				mycount = float(info[myindex]) / mylen
				mystr = mystr + "\t" + str(mycount)
				myindex = myindex + 1
			open_out.write(mystr + "\n")
		# if exist cluster length info
		else:
			# debug
			print("No length info!\t" + myclust)
	# foreach line
	open_file.close()
	open_out.close()
# normalization


#==============================================================
###########  Main processing ############
#==============================================================
def main():

	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start abundance_RPK.py -i " + values.i + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	cluster = collect_cluster_info (values.l, values.t)
	sys.stderr.write("Get cluster info ......done\n")
	
	### normalization ###
	sys.stderr.write("Normalize abundance......starting\n")
	normalization (cluster, values.i, values.o)
	sys.stderr.write("Normalizae abundance......done\n")


	sys.stderr.write("### Finish abundance_RPK.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

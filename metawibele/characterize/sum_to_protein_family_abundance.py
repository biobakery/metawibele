#!/usr/bin/env python

"""
MetaWIBELE: sum_to_protein_family_abundance module
Sum over protein families based on reads mapped to gene catalog of nucleotide representatives

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
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Sum over protein families based on reads mapped to gene catalog of nucleotide representatives
"""

def get_args ():	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', help='input gene catalog abundance', required=True)
	parser.add_argument('-t', help='specify the type of data tables, e.g. count | RPK | relab | CPM', default="count")
	parser.add_argument('-o', help='output abundance file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_protein_cluster_info (clust_file):	# discovery_cohort.proteins.clust
	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	myrep = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search("cluster=([\d]+)", line)
			#myclust = "Cluster_" + mym.group(1)
			myclust = mym.group(1)
			mym = re.search("^>([^\;]+)", line)
			myrep = mym.group(1)
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myclust
	# foreach line
	open_file.close()
	return cluster
# function collect_gene_cluster_info



#==============================================================
# sum up to cluster abundance
#==============================================================
def assign_counts (pep_cluster, data_type, infile, outfile):
	counts = {}
	titles = {}
	samples = {}
	outs = {}
	open_file = open(infile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^ID", line):
			myindex = 1
			while myindex < len(info):
				item = info[myindex]
				titles[myindex] = item
				samples[item] = ""
				myindex = myindex + 1
			# foreach item
			continue
		# if title
		myid = info[0]
		if not myid in pep_cluster:
			continue
		myclust = pep_cluster[myid]
		if not myclust in outs:
			outs[myclust] = {}
		myindex = 1
		while myindex < len(info):
			mycount = info[myindex]
			mys = titles[myindex]
			if not mys in outs[myclust]:
				outs[myclust][mys] = 0
			if data_type == "count":
				outs[myclust][mys] = outs[myclust][mys] + int(mycount)
			else:
				outs[myclust][mys] = outs[myclust][mys] + float(mycount)
			myindex = myindex + 1
		# foreach sample
	# foreach line
	open_file.close()

#	for pepid in pep_cluster.keys(): # foreach protein family
#		for member in pep_cluster[pepid].keys():
#			myclust = pep_cluster[pepid][member]
#			if member in counts:
#				for mys in counts[member].keys():
#					if not myclust in outs:
#						outs[myclust] = {}
#					if not mys in outs[myclust]:
#						outs[myclust][mys] = 0
#					outs[myclust][mys] = outs[myclust][mys] + int(counts[member][mys])
#					#outs[myclust][mys][member] = int(counts[member][mys])
#				# foreach sample
#			# if geneid
#			#else:
#				# debug
#			#	print("No count info for this member\t" + member + "\t" + myclust)
#		# foreach member
#	# foreach protein cluster

	open_out = open(outfile, "w")
	mytitle = "ID\t" + "\t".join(sorted(samples.keys()))
	open_out.write(mytitle + "\n")
	for myclust in sorted(outs.keys()):
		mystr = "Cluster_" + myclust
		for mys in sorted(samples.keys()):
			mycount = 0
			if mys in outs[myclust]:
				mycount = outs[myclust][mys]
				#for myid in outs[myclust][mys]:
				#	mycount = mycount + outs[myclust][mys]
			mystr = mystr + "\t" + str(mycount)
		# foreach samples
		open_out.write(mystr + "\n")
	# foreach protein family
	open_out.close()
# assign_counts


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start sum_to_protein_family_abundance.py -i " + values.i + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	pep_cluster = collect_protein_cluster_info (config.protein_family)
	sys.stderr.write("Get cluster info ......done\n")

	### assign counts to protein families ###
	sys.stderr.write("\nAssign counts to protein families ......starting\n")
	assign_counts (pep_cluster, values.t, values.i, values.o)
	sys.stderr.write("\nAssign counts to protein families ......done\n")

	sys.stderr.write("### Finish sum_to_protein_family_abundance.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

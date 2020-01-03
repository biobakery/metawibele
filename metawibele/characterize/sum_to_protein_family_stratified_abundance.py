#!/usr/bin/env python3

"""
MetaWIBELE: sum_to_protein_family_stratified_abundance module
Sum over stratified protein families based on reads mapped to gene catalog of nucleotide representatives

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
	parser.add_argument('-m', help='input taxonomy file for gene catalogs', required=True)
	parser.add_argument('-t', help='threshold for taxon presence (e.g. >X% of genes non-zero in RNA)', default=0.5)
	parser.add_argument('-o', help='output stratified abundance file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_protein_cluster_info (clust_file):	
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
# collect taxa info
#==============================================================
def collect_taxa_info (taxa_file):
	taxa_levels = ["Species"]
	taxa = {}
	taxa_num = {}
	titles = {}
	open_file = open(taxa_file, "r")
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
		myid = info[0]
		mytaxa = "NA"
		mylevel = "NA"
		if "taxa_name" in titles:
			mytaxa = info[titles["taxa_name"]]
			if "taxa_rank" in titles:
				 mylevel = info[titles["taxa_rank"]]
		else:
			if "detail" in titles:
				mytaxa = info[titles["detail"]]
		if mylevel != "NA":
			if not mylevel in taxa_levels:
				continue
		if mytaxa == "NA" or mytaxa == "Unclassified":
			continue
		taxa[myid] = mytaxa
		if not mytaxa in taxa_num:
			taxa_num[mytaxa] = {}
		taxa_num[mytaxa][myid] = ""
	# foreach line
	open_file.close()

	return taxa, taxa_num
# collect_taxa_info


#==============================================================
# sum up to family abundance
#==============================================================
def sum_abundance (pep_cluster, taxa, taxa_num, flt_presence, infile, outfile):
	titles = {}
	samples = []
	taxa_presence = {}
	outs = {}
	open_file = open(infile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 1
	while myindex < len(info):
		item = info[myindex]
		titles[myindex] = item
		samples.append(item)
		myindex = myindex + 1
	# foreach item
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if not myid in pep_cluster:
			continue
		if myid in taxa:
			mytaxa = taxa[myid]
		else:
			continue
		myclust = pep_cluster[myid] + "|" + mytaxa
		if not mytaxa in taxa_presence:
			taxa_presence[mytaxa] = {}
		taxa_presence[mytaxa][myid] = ""
		if not myclust in outs:
			outs[myclust] = {}
		myindex = 1
		while myindex < len(info):
			mycount = info[myindex]
			mys = titles[myindex]
			if not mys in outs[myclust]:
				outs[myclust][mys] = 0
			if re.search("\.", mycount):
				outs[myclust][mys] = outs[myclust][mys] + float(mycount)
			else:
				outs[myclust][mys] = outs[myclust][mys] + int(mycount)
			myindex = myindex + 1
		# foreach sample
	# foreach line
	open_file.close()

	# threshold for species presence (e.g. >X% of genes non-zero in RNA)
	refined_taxa = {}
	for mytaxa in taxa_presence:
		mynum = len(taxa_presence[mytaxa].keys())
		if mytaxa in taxa_num:
			mytotal = len(taxa_num[mytaxa].keys())
			myper = float(mynum) / float(mytotal)
			if myper >= float(flt_presence):
				refined_taxa[mytaxa] = ""

	# output abundance
	samples = list(set(samples)) 
	open_out = open(outfile, "w")
	mytitle = "ID\t" + "\t".join(samples)
	open_out.write(mytitle + "\n")
	for myclust in sorted(outs.keys()):
		mystr = "Cluster_" + myclust
		tmp = myclust.split("|")
		if not tmp[1] in refined_taxa:
			continue
		for mys in samples:
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


	sys.stderr.write("### Start sum_to_protein_family_stratified_abundance.py -i " + values.i + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	pep_cluster = collect_protein_cluster_info (config.protein_family)
	taxa, taxa_num = collect_taxa_info (values.m)
	sys.stderr.write("Get cluster info ......done\n")

	### sum abundance to protein families ###
	sys.stderr.write("\nAssign counts to protein families ......starting\n")
	sum_abundance (pep_cluster, taxa, taxa_num, values.t, values.i, values.o)
	sys.stderr.write("\nAssign counts to protein families ......done\n")

	sys.stderr.write("### Finish sum_to_protein_family_stratified_abundance.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

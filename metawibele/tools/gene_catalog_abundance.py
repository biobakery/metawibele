#!/usr/bin/env python

"""
MetaWIBELE: gene_catalog_abundance module
Summary abundance table of gene catalog

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
Summary abundance table of gene catalog
"""
	
def get_args (): 	
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', help='input the working directory path of counts files', required=True)
	parser.add_argument('-s', help='specify the suffix of counts file, e.g. sort.bed', default="sort.bed")
	parser.add_argument('-c', help='input non-redundant gene catalog clustering info file', required=True)
	parser.add_argument('-o', help='output cluster distribution file', required=True)
	values = parser.parse_args()
	return values
# get_args	


#==============================================================
# collect cluster info
#==============================================================
def collect_gene_cluster_info (clust_file):	# discovery_cohort.genes.clust
	cluster = {}
	open_file = open(clust_file, "r")
	myrep = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			#mym = re.search("cluster=([\d]+)", line)
			#myclust = mym.group(1)
			mym = re.search("^>([^\;]+)", line)
			myrep = mym.group(1)
			cluster[myrep] = myrep
			continue
	# foreach line
	open_file.close()
	return cluster
# function collect_gene_cluster_info


#==============================================================
# assign counts to peptide families info
#==============================================================
def collect_counts(map_path, extension, gene_cluster):
	counts = {}
	mysample = {}
	
	'''
	samples = {}
	open_file = open(sample_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		samples[line.split("\t")[0]] = ""
	# foreach sample
	open_file.close()
	'''

	filelist = utilities.find_files(map_path, extension, None)
	for myfile in filelist:
		mym = re.search("([^\/]+)$", myfile)
		sample = mym.group(1)
		sample = re.sub("." + extension, "", sample)
		if not os.path.isfile(myfile):
			config.logger.info ("ERROR! File not exist: " + myfile)
			continue
		mysample[sample] = ""
		open_file = open(myfile, "r")
		for line in open_file:
			line = line.strip()
			if not len(line):
				continue
			if re.search("^#", line):
				continue
			if re.search("^Geneid", line):
				continue
			info = line.split("\t")
			myid = info[0]
			if not myid in gene_cluster:	# not specified genes
				continue
			mycount = info[-1]
			if mycount == str(0):	# no counts
				continue
			if not myid in counts:
				counts[myid] = {}
			counts[myid][sample] = mycount
		# foreach line
		open_file.close()
	# foreach samplelist
	
	return counts, mysample
# collect_counts


# assign_counts
def assign_counts (counts, samples, outfile):
	open_out = open(outfile, "w")
	mytitle = "ID\t" + "\t".join(sorted(samples.keys())) 
	open_out.write(mytitle + "\n")
	for myclust in sorted(counts.keys()):
		mystr = myclust
		for mys in sorted(samples.keys()):
			if mys in counts[myclust].keys():
				mycount = counts[myclust][mys]
			else:
				mycount = "0"
			mystr = mystr + "\t" + str(mycount)
		# foreach sample
		open_out.write(mystr + "\n")
	# foreach peptide family
	open_out.close()
# assign_counts


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	### get arguments ###
	values = get_args()

	config.logger.info ("#### Start gene_catalog_abundance step ####")


	### collect cluster info ###
	config.logger.info ("Get cluster info ......starting")
	gene_cluster = collect_gene_cluster_info (values.c)
	config.logger.info ("Get cluster info ......done")
	
	### collect counts info ###
	config.logger.info ("Get counts info ......starting")
	counts, samples = collect_counts(values.p, values.s, gene_cluster)
	config.logger.info ("Get counts info ......done")

	### assign counts to peptide families ###
	config.logger.info ("Assign counts to peptide families ......starting")
	assign_counts (counts, samples, values.o)
	config.logger.info ("Assign counts to peptide families ......done")

	config.logger.info ("### Finish gene_catalog_abundance step ####")

# end: main

if __name__ == '__main__':
	main()

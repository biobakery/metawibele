#!/usr/bin/env python
"""
MetaWIBELE: abundance_annotator module
Summary the abundance for each taxon of each protein family across samples

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
import math

try:
	from metawibele import utilities
	from metawibele import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Summary the abundance for proteinfamilies/protein 
"""


def get_args():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--abundance",
	                    help='input family abundance file',
	                    required=True)
	parser.add_argument('-t', "--type",
	                    help='specify the type of abundance table, e.g. DAN_abundance | RNA_abundance | RNA-ratio_abundance>',
	                    required=True)
	parser.add_argument('-f', "--flag",
	                    help='specify the type of features',
	                    choices=["protein", "protein_family"],
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output annotated file',
	                    required=True)
	values = parser.parse_args()

	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_cluster_info (clust_file, flag):  # discovery_cohort.peptides.clust
	cluster = {}
	cluster_mem = {}
	open_file = open(clust_file, "r")
	myclust = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			mym = re.search("cluster=([\d]+)", line)
			myclust_id = "Cluster_" + mym.group(1)
			cluster[myclust_id] = myclust
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster_mem[myid] = myclust + "\t" + myclust_id
    # foreach line
	open_file.close()

	if flag == "protein":
		return cluster_mem
	if flag == "protein_family":
		return cluster
# function collect_cluster_info


# ==============================================================
# collect cluster abundance
# ==============================================================
def collect_cluster_abundance (abufile, sample_info):
	abundance = {}
	sample_num = {}
	titles = {}
	open_file = open(abufile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[info.index(item)] = item
	sample_num["combined"] = {}
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myindex = 1
		while myindex < len(info):
			myvalue = info[myindex]
			mys = titles[myindex]
			sample_num["combined"][mys] = ""
			mytype = "NA"
			if mys in sample_info:
				mytype = sample_info[mys]
			if mytype != "NA":
				if not mytype in sample_num:
					sample_num[mytype] = {}
				sample_num[mytype][mys] = ""
			if myvalue != "NA" and myvalue != "NaN" and myvalue != "nan" and myvalue != "inf" and myvalue != "Inf":
				myvalue = float(myvalue)
				if myvalue > float(config.abundance_detection_level):
					if not myid in abundance:
						abundance[myid] = {}
					if not "combined" in abundance[myid]:
						abundance[myid]["combined"] = []
					abundance[myid]["combined"].append(myvalue)
					if mytype != "NA":
						if not mytype in abundance[myid]:
							abundance[myid][mytype] = []
						abundance[myid][mytype].append(myvalue)
			myindex = myindex + 1
	# foreach line
	open_file.close()
	return abundance, sample_num
# collect_cluster_abundance


# ==============================================================
# output info
# ==============================================================
def output_info (cluster, sample_num, abundance, myflag, flag, outfile):
	myID = utilities.PROTEIN_FAMILY_ID
	if flag == "protein":
		myID = utilities.PROTEIN_ID
	open_out = open(outfile, "w")
	open_out.write(myID + "\ttype\tdetail\tmean_abundance\tmean_prevalent_abundance\tprevalence\ttotal_abundance\ttotal_samples\n")
	for myid in sorted(cluster.keys()):
		if myid in abundance:
			for mytype in sorted(abundance[myid].keys()):
				mynum = len(abundance[myid][mytype])
				mytotal = sum(abundance[myid][mytype])
				mypre = float(mynum) / float(len(sample_num[mytype].keys()))
				mymean_beta = utilities.mean(abundance[myid][mytype])
				mymean = float(mytotal) / float(len(sample_num[mytype].keys()))
				myflag1 = re.sub("_abundance", "-" + mytype + "_abundance", myflag)
				myflag1 = re.sub("-combined_abundance", "_abundance", myflag1)
				myflag1 = re.sub("non_", "non-", myflag1)
				mystr = myid + "\t" + myflag1 + "\t" + str(mymean_beta) + "\t" + str(mymean) + "\t" + str(mymean_beta) + "\t" + str(mypre) + "\t" + str(mytotal) + "\t" + str(len(sample_num[mytype].keys()))
				open_out.write(mystr + "\n")
				
				myflag2 = re.sub("abundance", "prevalence", myflag1)
				mystr = myid + "\t" + myflag2 + "\t" + str(mypre) + "\t" + str(mymean) + "\t" + str(mymean_beta) + "\t" + str(mypre) + "\t" + str(mytotal) + "\t" + str(len(sample_num[mytype].keys()))
				open_out.write(mystr + "\n")
		else:
			mystr = myid + "\t" + myflag + "\t" + "NA" + "\t" + "NaN" + "\t" + "NaN" + "\t" + "NaN" + "\t" + "NaN" + "\t" + str(len(sample_num["combined"].keys()))
			open_out.write(mystr + "\n")
			myflag1 = re.sub("abundance", "prevalence", myflag)
			mystr = myid + "\t" + myflag1 + "\t" + "NA" + "\t" + "NaN" + "\t" + "NaN" + "\t" + "NaN" + "\t" + "NaN" + "\t" + str(len(sample_num["combined"].keys()))
			open_out.write(mystr + "\n")
	open_out.close()
# combine info


# ==============================================================
###########  Main processing ############
# ==============================================================
if __name__ == '__main__':
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start abundance_annotator.py -a " + values.abundance + " ####\n")

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	sample_info = utilities.sample_info (config.metadata, config.study)
	cluster = collect_cluster_info (config.protein_family, values.flag)
	abundance, sample_num = collect_cluster_abundance (values.abundance, sample_info)
	sys.stderr.write("Get info ......done\n")

	### convert info ###
	sys.stderr.write("Convert info ......starting\n")
	output_info (cluster, sample_num, abundance, values.type, values.flag, values.output)
	sys.stderr.write("Output info ......done\n")

	sys.stderr.write("### Finish abundance_annotator.py ####\n\n\n")

# end: main
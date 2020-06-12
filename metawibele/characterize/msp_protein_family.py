#!/usr/bin/env python

"""
MetaWIBELE: msp_protein_family module
Summary taxonomy novelty based on MSP for each protein family

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
Summary taxonomy novelty based on MSP for each protein family
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-m', "--msp",
	                    help='input MSP taxonomy annotation file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output taxonomy novelty annotation detailed file',
	                    required=True)
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
	myclust_id = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			mym = re.search("cluster=([\d]+)", line)
			myclust_id = "Cluster_" + mym.group(1)
			if not myclust in cluster:
				cluster[myclust] = {}
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myclust][myid] = myclust_id
	# foreach line
	open_file.close()
	return cluster
# function collect_protein_cluster_info

def collect_gene_cluster_info (clust_file):	
	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myclust
	# foreach line
	open_file.close()
	return cluster
# function collect_gene_cluster_info


#==============================================================
# assign annotation to protein families info
#==============================================================
def collect_annotation(mspfile):	# HMP2_MSPminer_annotation.taxonomy.tsv
	annotation = {}
	titles = {}
	open_file = open(mspfile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	ann_title = "\t".join(info[1:len(info)])
	myindex = 0
	while myindex < len(info):
		item = info[myindex]
		titles[item] = myindex
		myindex = myindex + 1
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if myid in annotation:
			continue
		mytype = "Novel_taxonomy"
		myinfo = "\t".join(info[1:len(info)])
		taxa = info[titles["taxa_rank"]]
		source = info[titles["source"]]
		if re.search("UniRef90_unclassified", source):
			if taxa == "Unclassified":
				myclass = "unclassified"
			else:
				myclass = "novel_classified"
		else:
			if taxa == "Unclassified":
				# debug
				print("UnRef90\tunclassified\t" + myid)
				myclass = "unclassified"
			else:
				myclass = "known_classified"
		if not myid in annotation:
			annotation[myid] = []
		annotation[myid].append(mytype + "\t" + myclass + "\n" + myinfo + "\n" + taxa)
	# foreach line
	open_file.close()
	return annotation, ann_title
# collect_annotation

# assign_annotation
def assign_annotation (propt_cluster, annotation, ann_title, outfile_detail):
	outs_info = {}
	outs_ORF = {}
	titles = {}
	ann_info = {}
	cluster = {}
	types = ["novel_classified", "known_classified", "unclassified"]
	tmp = ann_title.split("\t")
	for item in tmp:
		titles[item] = tmp.index(item)
	for pepid in sorted(propt_cluster.keys()): # foreach protein family
		flag = 0
		detail_info = {}
		mytotal = 0
		clust_id = ""
		for member in propt_cluster[pepid].keys():
			gene_id = member
			clust_id = propt_cluster[pepid][member]
			if not gene_id in annotation:	# no corresponding annotation
				# debug
				print("No MSP-based taxonomy annotation info for the member!\t" + gene_id)
				continue
			for myall in annotation[gene_id]:
				myid, myinfo, mynote = myall.split("\n")
				mytype, ann = myid.split("\t")
				tmp_info = myinfo.split("\t")
				if not clust_id in ann_info:
					ann_info[clust_id] = {}
				if not member in ann_info[clust_id]:
					ann_info[clust_id][member] = ann + "\t" + mynote
				if not ann in detail_info:
					detail_info[ann] = {}
				detail_info[ann][member] = ""
			# foreach type
		# foreach member
		cluster[clust_id] = ""

		# percentage
		detail = {}
		detail_ORF = {}
		mytotal = len(propt_cluster[pepid].keys())
		for mytype in detail_info.keys():
			mynum = len(detail_info[mytype].keys())
			detail[mytype] = mynum
			for member in detail_info[mytype].keys():
				detail_ORF[member] = mytype
		# foreach type
		
		mystr = "Novel_taxonomy"
		for mytype in types:
			mynum = 0
			if mytype in detail:
				mynum = detail[mytype]
			if mytype == "novel_classified":
				myper = float(mynum) / float(mytotal)
				mystr = mystr + "\t" + str(myper)
			mystr = mystr + "\t" + str(mynum)
		# foreach type
		outs_info[clust_id] = mystr + "\t" + str(mytotal)
		
		for member in detail_ORF.keys():
			mytype = detail_ORF[member]
			mystr1 = "Novel_taxonomy"
			for mytype1 in types:
				mynum1 = 0
				if mytype1 == mytype:
					mynum1 = 1
				if mytype1 == "novel_classified":
					myper1 = float(mynum1)
					mystr1 = mystr1 + "\t" + str(myper1)
				mystr1 = mystr1 + "\t" + str(mynum1)
			# foreach type
			outs_ORF[member] = mystr1 + "\t" + str(1)
		# foreach member
	# foreach protein cluster
	
	#### output info ####
	outfile1 = re.sub(".detail.tsv", ".ORF.detail.tsv", outfile_detail)
	open_out = open(outfile_detail, "w")
	open_out1 = open(outfile1, "w")
	open_out.write(utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\t" + "novel_classified\tknown_classified\tunclassified\ttotal_member" + "\n")
	open_out1.write(utilities.PROTEIN_ID + "\ttype\tdetail\t" + "novel_classified\tknown_classified\tunclassified\ttotal_member" + "\n")
	for myclust in sorted(outs_info.keys()):
		open_out.write(myclust + "\t" + outs_info[myclust] + "\n")
	# foreach cluster
	open_out.close()
	for mypep in sorted(outs_ORF.keys()):
		open_out1.write(mypep + "\t" + outs_ORF[mypep] + "\n")
	# foreach protein family
	open_out1.close()

# assign_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start msp_protein_family.py -m " + values.msp + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	propt_cluster = collect_protein_cluster_info (config.protein_family)
	sys.stderr.write("Get cluster info ......done\n")
	
	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	annotation, ann_title = collect_annotation(values.msp)
	sys.stderr.write("Get annotation info ......done")

	### assign annotation to protein families ###
	sys.stderr.write("\nAssign annotation to protein families ......starting\n")
	assign_annotation (propt_cluster, annotation, ann_title, values.output)
	sys.stderr.write("\nAssign annotation to protein families ......done\n")

	sys.stderr.write("### Finish msp_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

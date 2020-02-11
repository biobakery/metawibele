#!/usr/bin/env python

"""
MetaWIBELE: extract_sequences module
Extract gene and protein sequence based on clusterID

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


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Extract gene sequence based on clusterID 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-c', help='input cluster info file', required=True)
	parser.add_argument('-p', help='input protein cluster info file', required=True)
	parser.add_argument('-g', help='input gene cluster info file', required=True)
	parser.add_argument('-s', help='input gene sequence file', required=True)
	parser.add_argument('-q', help='input protein sequence file', required=True)
	parser.add_argument('-o', help='output cluster sequence file', required=True)
	values = parser.parse_args()
	return values
# get_args


# remove duplicates in list
def remove_duplicate (duplicate): 
	final_list = [] 
	for num in duplicate: 
		if num not in final_list: 
			final_list.append(num) 
	return final_list 
# remove_duplicate


#==============================================================
# collect cluster info
#==============================================================
def collect_pep_cluster_info (clust_file):	
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
			if not myclust_id in cluster:
				cluster[myclust_id] = myclust
			continue
	# foreach line
	open_file.close()
	return cluster
# function collect_pep_cluster_info

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
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myclust
	# foreach line
	open_file.close()
	return cluster
# function collect_gene_cluster_info


#==============================================================
# collect sequence info 
#==============================================================
def collect_seq_info (seq_file):
	seqs = {}
	open_file = open(seq_file, "r")
	myid = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([\S]+)", line)
			myid = mym.group(1)
			if not myid in seqs:
				seqs[myid] = ""
			continue
		seqs[myid] = seqs[myid] + line
	# foreach line
	open_file.close()
	return seqs
# collect_seqs_info


#==============================================================
# collect cluster ID info
#==============================================================
def collect_ID_info (cluster_file):
	cluster = []
	titles = {}
	flags = {}
	open_file = open(cluster_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 0
	for item in info:
		titles[item] = info.index(item)
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[titles["familyID"]]
		myname = myid
		if myid in flags:
			continue
		flags[myid] = ""
		if "unirefID" in titles:
			myuniref = info[titles["unirefID"]]
			if "map_type" in titles:
				mymap = info[titles["map_type"]]
				if re.search("weak_homology", mymap) or re.search("worse_homology", mymap):
					myuniref = "NA"
			if myuniref != "NA":
				myname = myid + "|" + myuniref
		cluster.append(myid + "\t" + myname) 
	# foreach line
	open_file.close()
	return cluster
# collect_ID_info


#==============================================================
# assign sequence to cluster
#==============================================================
def assign_gene_seq (clusterID, pep_cluster, gene_cluster, seqs, outfile):
	outfile1 = outfile + ".gene.fna"
	open_file = open(outfile1, "w")
	for mytmp in clusterID: 
		myclust, myname = mytmp.split("\t")
		if myclust in pep_cluster:
			myid = pep_cluster[myclust]	 # protein id
			if not myid in gene_cluster:  # no corresponding gene cluster
				print("Peptide ID has no corresponding gene cluster!\t" + myid)
				continue
			gene_id = gene_cluster[myid]
			if gene_id in seqs:
				open_file.write(">" + myname + "\n" + seqs[gene_id] + "\n")
		# cluster name
	# foreach peptide cluster with DA
	open_file.close()	
# assign_gene_seq

def assign_protein_seq (clusterID, pep_cluster, seqs, outfile):
	outfile1 = outfile + ".protein.faa"
	open_file = open(outfile1, "w")
	for mytmp in clusterID:
		myclust, myname = mytmp.split("\t")
		if myclust in pep_cluster:
			myid = pep_cluster[myclust]	 # protein id
			if myid in seqs:
				open_file.write(">" + myname + "\n" + seqs[myid] + "\n")
		# cluster name
	# foreach peptide cluster with DA
	open_file.close()	
# assign_protein_seq



#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start extract_sequences.py -p " + values.p + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	sys.stderr.write("Collect total cluster info ......\n")
	pep_cluster = collect_pep_cluster_info (values.p)
	gene_cluster = collect_gene_cluster_info (values.g)
	sys.stderr.write("Collect annotated cluster info ......\n")
	cluster = collect_ID_info (values.c)
	sys.stderr.write("Collect sequence info ......\n")
	gene_seqs = collect_seq_info (values.s)
	pep_seqs = collect_seq_info (values.q)
	sys.stderr.write("Get info ......done\n")

	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign sequence to peptide families ......starting\n")
	assign_gene_seq (cluster, pep_cluster, gene_cluster, gene_seqs, values.o)
	assign_protein_seq (cluster, pep_cluster, pep_seqs, values.o)
	sys.stderr.write("\nAssign sequence to peptide families ......done\n")

	sys.stderr.write("### Finish extract_sequences.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

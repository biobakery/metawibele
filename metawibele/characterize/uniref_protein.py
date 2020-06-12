#!/usr/bin/env python

"""
MetaWIBELE: uniref_protein module
Extract proteins annotated by UniRef DBs

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
Extract predicted peptides annotated by UniRef DBs
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-m', "--mapping",
	                    help='input the annotated stat file based on UniRef DB',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output the annotation results',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_peptide_cluster_info (clust_file):  # discovery_cohort.peptides.clust
	member = {}
	open_file = open(clust_file, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		member[myid] = ""
    # foreach line
	open_file.close()
	return member
# function collect_peptide_cluster_info


#==============================================================
# collect Pfam info
#==============================================================
def collect_pfam_info (pfamfile):	# Pfam_ann.tsv
	pfam = {}
	open_file = open(pfamfile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^Pfam", line):
			continue
		pfam[info[0]] = info[1]
	# foreach line
	open_file.close()
	return pfam
# collect_pfam_info


#==============================================================
# collect UniRef info
#==============================================================
def collect_uniref_info (uniref, pfam, hits):
	uniref_info = {}
	titles = {}
	#items = ["UniRefID", "Tax", "TaxID", "Description", "Organism", "GO(MF)", "GO(CC)", "Subcellular_location", "Transmembrane", "Signal_peptide", "Pfam"]
	open_file = open(uniref, "r")
	line = open_file.readline()
	line = line.strip()
	uniref_title = line
	uniref_title = re.sub("^ID", "UniRefID", uniref_title)
	uniref_title = uniref_title + "\tPfam_desc"
	info = line.split("\t")
	myindex = 0
	while myindex < len(info):
		item = info[myindex]
		titles[item] = myindex
		myindex = myindex + 1
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if not myid in hits:
			continue
		mystr = line
		mystr = re.sub("root\tNA", "root\t1", mystr)
		mypfam = info[titles["Pfam"]]
		if mypfam == "NA":
			mystr = mystr + "\tNA"
		else:
			tmp = mypfam.split(";")
			mypfam_info = ""
			for item in tmp:
				if item in pfam:
					mypfam_info = mypfam_info + pfam[item] + ";"
				else:
					mypfam_info = mypfam_info + "NA;"
			# foreach pfam
			mypfam_info = re.sub(";$", "", mypfam_info)
			mystr = mystr + "\t" + mypfam_info
		if not myid in uniref_info:
			uniref_info[myid] = []
			# only use the annotation of representatives
			uniref_info[myid].append(mystr)
	# foreach line
	open_file.close()
	return uniref_info, uniref_title
# collect_uniref_info


#==============================================================
# collect UniRef mapping info
#==============================================================
def collect_uniref_mapping (mapfile, member):	# PRISM_peptides.clust.rep.uniref90.stat.all.tsv
	titles = {}
	mapping = {}
	hits = {}
	open_file = open(mapfile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myid = info[titles[utilities.PROTEIN_ID]]
		if not myid in member:
			# debug
			#print("Not members of specified clusters!\t" + myid)
			continue
		myuniref = info[titles["subject"]]
		hits[myuniref] = ""
		query_type = info[titles["query_type"]]
		mutual_type = info[titles["mutual_type"]]
		mapping[myid] = myuniref + "\t" + query_type + "\t" + mutual_type
	# foreach line
	open_file.close()

	return mapping, hits
# collect_uniref_mapping


#==============================================================
# extract ann info
#==============================================================
def extract_annotation_info (uniref_info, uniref_title, mapping, outfile):
	open_out = open(outfile, "w")
	open_out.write("name\t" + uniref_title + "\n")
	item_num = len(uniref_title.split("\t"))
	for myid in sorted(mapping.keys()):
		info = mapping[myid].split("\t")
		myuniref = info[0]
		query_type = info[1]
		mutual_type = info[2]
		uniref = "UniRef90_unknown"
		if mutual_type == "high_confidence":
			uniref = myuniref
		if uniref in uniref_info:
			for item in uniref_info[uniref]:
				open_out.write(myid + "\t" + item + "\n")
		else:
			mystr = myid + "\t" + uniref
			mynum = 1
			while mynum < item_num:
				mystr = mystr + "\tNA"
				mynum = mynum + 1
			open_out.write(mystr + "\n")
	# foreach line
	open_out.close()
# function extract_annotation_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start uniref_protein.py -m " + values.mapping+ " ####\n")
	
	### collect uniref and annotation info ###
	sys.stderr.write("Get UniRef DB and annotation info ......starting\n")
	member = collect_peptide_cluster_info (config.protein_family)
	mapping, hits = collect_uniref_mapping (values.mapping, member)
	pfam = collect_pfam_info (config.pfam_database)
	uniref_info, uniref_title = collect_uniref_info (config.uniref_database, pfam, hits)
	extract_annotation_info (uniref_info, uniref_title, mapping, values.output)
	sys.stderr.write("Get UniRef DB and annotation info ......done\n")

	sys.stderr.write("### Finish uniref_protein.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

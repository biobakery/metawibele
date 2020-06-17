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
	from metawibele.common import utils
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
	for line in utils.gzip_bzip2_biom_open_readlines (pfamfile): 
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^Pfam", line):
			continue
		pfam[info[0]] = info[1]
	# foreach line
	
	return pfam
# collect_pfam_info


#==============================================================
# collect UniRef info
#==============================================================
def collect_basic_info (uniref_list, hits):
	uniref_info = {}
	uniref_info_tmp = {}
	flags = {}
	names = [] # ["UniRefID", "Protein_names", "Gene_names", "UniProtKB", "Tax", "TaxID", "Rep_Tax", "Rep_TaxID", "GO", "KO", "eggNOG", "Pfam", "Level4EC"]
	hits_items = set(sorted(hits.keys()))
	for myfile in uniref_list:
		myfile = myfile.strip()
		if not len(myfile):
			continue
		if re.search("^#", myfile):
			continue
		if not os.path.isfile(myfile):
			print("File does not exist! " + myfile)
			continue
		myname = os.path.basename(myfile)
		if not re.search("^map_", myname):
			continue
		if not re.search("uniref90", myname) and not re.search("uniref50", myname):
			continue
		myname = re.sub("map_", "", myname)
		if re.search("^[\S]+_uniref", myname):
			myname = re.sub("_uniref[\S]+$", "", myname)
			if myname == "go":
				myname = "GO"
			if myname == "ko":
				myname = "KO"
			if myname == "eggnog":
				myname = "eggNOG"
			if myname == "pfam":
				myname = "Pfam"
			if myname == "level4ec":
				myname = "Level4EC"
		if re.search("^uniref[\d]+_name", myname):
			myname = "Protein_names"
		if not myname in flags:
			flags[myname] = ""
			names.append(myname)
		else:
			continue
		for line in utils.gzip_bzip2_biom_open_readlines (myfile): 
			line = line.strip()
			if not len(line):
				continue
			if re.search("^#", line) :
				continue
			info = line.split("\t")
			myvalue = info[0]
			myset = set(sorted(info[1:len(info)]))
			myoverlap = myset.intersection(hits_items)
			for mykey in myoverlap:
				if not mykey in uniref_info_tmp:
					uniref_info_tmp[mykey] = {}
				if not myname in uniref_info_tmp[mykey]:
					uniref_info_tmp[mykey][myname] = myvalue
				else:
					uniref_info_tmp[mykey][myname] = uniref_info_tmp[mykey][myname] + ";" + myvalue

			#myindex = 1
			#while myindex < len(info):
			#	mykey = info[myindex]
			#	if not mykey in hits:
			#		myindex = myindex + 1
			#		continue
			#	if not mykey in uniref_info_tmp:
			#		uniref_info_tmp[mykey] = {}
			#	if not myname in uniref_info_tmp[mykey]:
			#		uniref_info_tmp[mykey][myname] = myvalue
			#	else:
			#		uniref_info_tmp[mykey][myname] = uniref_info_tmp[mykey][myname] + ";" + myvalue
			#	myindex = myindex + 1
		
		# foreach line
	# foreach dataset

	# collecting all info
	for myid in uniref_info_tmp.keys():
		mystr = myid
		for myname in names:
			if myname in uniref_info_tmp[myid]:
				mystr = mystr + "\t" + uniref_info_tmp[myid][myname]
			else:
				mystr = mystr + "\tNA"
		uniref_info[myid] = mystr
	uniref_info_tmp = {}

	return uniref_info, names
# function collect_expression_info


def collect_uniref_info (uniref, names, pfam):
	uniref_info = {}
	titles = {}
	uniref_title = "UniRefID\t" + "\t".join(names) + "\tPfam_desc"
	titles["UniRefID"] = 0
	myindex = 0
	while myindex < len(names):
		item = names[myindex]
		titles[item] = myindex + 1
		myindex = myindex + 1
	for myid in sorted(uniref.keys()):
		line = uniref[myid]
		info = line.split("\t")
		myid = info[0]
		mystr = line
		mystr = re.sub("root\tNA", "root\t1", mystr)
		if "Pfam" in titles:
			mypfam = info[titles["Pfam"]]
		else:
			mypfam = "NA"
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
	
	return uniref_info, uniref_title
# collect_uniref_info


#==============================================================
# collect UniRef mapping info
#==============================================================
def collect_uniref_mapping (mapfile, member):
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
		query_type = info[titles["query_type"]]
		mutual_type = info[titles["mutual_type"]]
		mapping[myid] = myuniref + "\t" + query_type + "\t" + mutual_type
		hits[myuniref] = ""
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
	uniref, names = collect_basic_info(config.uniref_database, hits)
	uniref_info, uniref_title = collect_uniref_info (uniref, names, pfam)
	extract_annotation_info (uniref_info, uniref_title, mapping, values.output)
	sys.stderr.write("Get UniRef DB and annotation info ......done\n")

	sys.stderr.write("### Finish uniref_protein.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

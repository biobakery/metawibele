#!/usr/bin/env python

"""
MetaWIBELE: mspminer_msp_uniref_annotation module
Assign UniRef annotation to MSPs 

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
Assign UniRef annotation to MSPs 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--annotation",
	                    help='input the uniref annotation file',
	                    required=True)
	parser.add_argument('-i', "--msp",
	                    help='input the MSPs summary information file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output summary file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# extract taxon for MSPs
#==============================================================
def collect_taxonomy_info (map_file):
	taxa_map = {}
	titles = {}
	open_file = open(map_file, "r")
	for line in utils.gzip_bzip2_biom_open_readlines (map_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^Taxon", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		# title line
		taxa_id = info[titles["Taxon"]]
		taxa_name = info[titles["Scientific_name"]]
		taxa_rank = info[titles["Rank"]]
		taxa_lineage = info[titles["Lineage"]]
		tmp = taxa_lineage.split("|")
		myid = taxa_name
		if re.search("__", tmp[-1]):
			mym = re.search("__([\S]+)", tmp[-1])
			myid = mym.group(1)
		taxa_map[myid] = taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage
    # foreach line
	
	return taxa_map
# collect_taxonomy_info


#==============================================================
# assign uniref annotation info to MSPs
#==============================================================
def assign_uniref_info (ann_file, msp_file, taxa_map, taxa_level, outfile):
	# collect info
	anns = {}
	ann_title = ""
	titles = {}
	open_file = open(ann_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search(utilities.PROTEIN_FAMILY_ID, line) or re.search(utilities.PROTEIN_ID, line):
			ann_title = line
			for item in info:
				titles[item] = info.index(item)
			continue
		myid = info[0]
		taxa_line = info[titles["taxa_lineage"]]
		if taxa_level != "no":
			mylevel = taxa_level
			if taxa_level == "Species":
				mylevel = "s__"
			if taxa_level == "Genus":
				mylevel = "g__"
			if taxa_level == "Family":
				mylevel = "f__"
			if taxa_level == "Order":
				mylevel = "o__"
			if taxa_level == "Class":
				mylevel = "c__"
			if taxa_level == "Phylum":
				mylevel = "p__"
			if taxa_level == "Kingdom":
				mylevel = "k__"
			if re.search(mylevel, taxa_line):
				mym = re.search("s__([^\|]+)", taxa_line)
				myname = mym.group(1)
				#myname = re.sub("_", " ", myname)
				if myname in taxa_map:
					taxa_id, taxa_name, taxa_rank, taxa_lineage = taxa_map[myname].split("\t")
					info[titles["taxa_id"]] = taxa_id
					info[titles["taxa_name"]] = taxa_name
					info[titles["taxa_rank"]] = taxa_rank
					info[titles["taxa_lineage"]] = taxa_lineage
				else:
					# debug
					print("No taxon info for this taxon name!\t" + myname)
		# if taxa_level
		line = "\t".join(info)
		anns[myid] = line
	# foreach line
	open_file.close()

	# assign info
	open_out = open(outfile, "w")
	open_file = open(msp_file, "r")
	titles = {}
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("msp_name\t", line):
			mystr = line + "\t" + ann_title
			open_out.write(mystr + "\n")
			for item in info:
				titles[item] = info.index(item)
			continue
		myid = info[titles["gene_name"]]
		if myid in anns:
			mystr = line + "\t" + anns[myid]
			open_out.write(mystr + "\n")
		else:
			# debug
			print("No corresponding unref annotation info!\t" + myid)
			#open_out.write(line + "\n")
	# foreach line
	open_out.close()
	open_file.close()
# func: assign_uniref_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start mspminer_msp_uniref_annotation.py -i " + values.msp + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	taxa_level = "no"
	taxa_map = collect_taxonomy_info (config.taxonomy_database)
	assign_uniref_info (values.annotation, values.msp, taxa_map, taxa_level, values.output)
	sys.stderr.write("Get info ......done\n")


	sys.stderr.write("### Finish mspminer_msp_uniref_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

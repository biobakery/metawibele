#!/usr/bin/env python3

"""
MetaWIBELE: antiSMASH_annotator module
Extract antiSMASH annotation based on UniRef ID

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
	parser.add_argument('-a', "--annotation",
	                    help='input UniRef annotation file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output the annotation results',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect BGC info
#==============================================================
def collect_bgc_info (bgcfile):
	bgc = {}
	bgc_ann = {}
	titles = {}
	open_file = open(bgcfile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^BGC_gene", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myuniref = info[titles["UniRef_ID"]]
		mybgc = info[titles["BGC_ID"]]
		mytype = info[titles["BGC_type"]]
		mydesc = info[titles["BGC_description"]]
		myclass = info[titles["BGC_class"]]
		mygene = info[titles["BGC_gene"]]
		bgc[mybgc] = mytype + "\t" + mydesc + "\t" + myclass
		if not myuniref in bgc_ann:
			bgc_ann[myuniref] = {}
		if not mybgc in bgc_ann[myuniref]:
			bgc_ann[myuniref] = {}
		if not mybgc in bgc_ann[myuniref]:
			bgc_ann[myuniref][mybgc] = {}
		bgc_ann[myuniref][mybgc][mygene] = ""
	# foreach line
	open_file.close()

	bgc_ann2 = {}
	for myid in bgc_ann.keys():
		if not myid in bgc_ann2:
			bgc_ann2[myid] = {}
		for mybgc in bgc_ann[myid].keys():
			mygene = ",".join(sorted(bgc_ann[myid][mybgc].keys()))
			bgc_ann2[myid][mybgc] = mygene

	return bgc, bgc_ann2

# func collect_bgc_info


#==============================================================
# collect UniRef annotation info
#==============================================================
def collect_uniref_annotation (annfile):
	titles = {}
	uniref = {}
	open_file = open(annfile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myflag = info[0]
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		map_type = info[titles["map_type"]]
		if map_type == "UniRef90_weak_homology" or map_type == "UniRef90_worse_homology":
			continue
		myuniref = info[titles["unirefID"]]
		if myuniref == "NA":
			continue
		uniref[myid] = myuniref
	# foreach line
	open_file.close()

	return uniref, myflag
# collect_uniref_annotation


#==============================================================
# get ann info
#==============================================================
def bgc_annotation (bgc, bgc_ann, uniref, id_flag, outfile):
	open_out = open(outfile, "w")
	open_out.write(id_flag + "\t" + "type\tdetail\tBGC_type\tBGC_description\tBGC_class\tBGC_gene" + "\n")
	for myid in sorted(uniref.keys()):
		myuniref = uniref[myid]
		if myuniref in bgc_ann:
			mydetial = "NA"
			mytype = "NA"
			mydesc = "NA"
			myclass = "NA"
			mygene = "NA"
			for mybgc in sorted(bgc_ann[myuniref].keys()):
				if mydetial == "NA":
					mydetial = mybgc
				else:
					mydetial = mydetial + ";" + mybgc
				if mygene == "NA":
					mygene = bgc_ann[myuniref][mybgc]
				else:
					mygene = mygene + ";" + bgc_ann[myuniref][mybgc]
				if mybgc in bgc:
					tmp = bgc[mybgc].split("\t")
					if mytype == "NA":
						mytype = tmp[0]
					else:
						mytype = mytype + ";" + tmp[0]
					if mydesc == "NA":
						mydesc = tmp[1]
					else:
						mydesc = mydesc + ";" + tmp[1]
					if myclass == "NA":
						myclass = tmp[2]
					else:
						myclass = myclass + ";" + tmp[2]
			# foreach BGC
			open_out.write(myid + "\t" + "antiSMASH_BGC" + "\t" + mydetial + "\t" + mytype + "\t" + mydesc + "\t" + myclass + "\t" + mygene + "\n")
	# foreach feature
	open_out.close()

# function bgc_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start antiSMASH_annotator.py -a " + values.annotation + " ####\n")
	
	### collect uniref and annotation info ###
	sys.stderr.write("Get antiSMASH DB and annotation info ......starting\n")

	bgc, bgc_ann = collect_bgc_info (config.antiSMASH_database)
	uniref, id_flag = collect_uniref_annotation (values.annotation)
	bgc_annotation (bgc, bgc_ann, uniref, id_flag, values.output)

	sys.stderr.write("Get antiSMASH DB and annotation info ......done\n")

	sys.stderr.write("### Finish antiSMASH_annotator.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

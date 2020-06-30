#!/usr/bin/env python

"""
MetaWIBELE: ddi_DOMINE_ann module
Assign domian description involved in DDI

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
Assign domian description involved in DDI
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--ddi-annotation",
	                    help='input DDIs summary file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output annotation summary file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect info
#==============================================================
def collect_DDI_info (DDI_file):	# summary_DOMINE_peptide.tsv
	DDIs = {}
	titles = {}
	open_file = open(DDI_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 0
	level = ["HC", "NA"]
	myflag = info[0]
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
		interact = info[titles["Interaction"]]
		pfam1 = info[titles["Pfam1_ID"]]
		pfam2 = info[titles["Pfam2_ID"]]
		ddi = pfam1 + "\t" + pfam2
		if not interact in level:
			continue
		if not myid in DDIs:
			DDIs[myid] = {}
		DDIs[myid][ddi] = ""
	# foreach line
	open_file.close()
	return myflag, DDIs
# function collect_DDI_info


def collect_pfam_info (pfamfile):	# Pfam_ann.txt 
	pfams = {}
	for line in utils.gzip_bzip2_biom_open_readlines (pfamfile):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		info = line.split("\t")
		pfam = info[0]
		ann = info[1]
		if not pfam in pfams:
			pfams[pfam] = {}
		pfams[pfam][ann] = ""
	# foreach line
	
	return pfams
#collect_pfam_info


def collect_pfam2go_info (pfam2go_file):	# Pfam2GO.txt 
	pfam2go = {}
	for line in utils.gzip_bzip2_biom_open_readlines (pfam2go_file): 
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		info = line.split("\t")
		pfam = info[0]
		go = info[1]
		ann = info[2]
		category = info[3]
		if category == "function":
			category = "MF"
		if category == "process":
			category = "BP"
		if category == "component":
			category = "CC"
		myann = ann + "(" + category + ":" + go + ")"
		if not pfam in pfam2go:
			pfam2go[pfam] = {}
		pfam2go[pfam][myann] = ""	
	# foreach line
	
	return pfam2go
# function collect_pfam2go_info


#==============================================================
# assign annotation to DDIs
#==============================================================
def assign_annotation (myflag, DDIs, pfams, pfam2go, outfile):
	anns = {}
	open_file = open(outfile, "w")
	open_file.write(myflag + "\tPfam1_ID\tPfam1_ann\tPfam2_ID\tPfam2_ann\n")
	for myclust in sorted (DDIs.keys()):
		for item in sorted(DDIs[myclust].keys()):
			pfam1, pfam2 = item.split("\t")
			ann1 = "NA"
			ann2 = "NA"
			if pfam1 in pfams:
				ann1 = ";".join(sorted(pfams[pfam1].keys()))
			if pfam2 in pfams:
				ann2 = ";".join(sorted(pfams[pfam2].keys()))
			open_file.write(myclust +"\t" + pfam1 + "\t" + ann1 + "\t" + pfam2 + "\t" + ann2 + "\n")
			myid = pfam1 + "\t" + ann1 + "\t" + pfam2 + "\t" + ann2

			if not myid in anns:
				anns[myid] = {}
			anns[myid][myclust] = ""
		# foreach type
	# foreach cluster
	open_file.close()
	
	outfile2 = re.sub(".tsv", ".GO.tsv", outfile)
	open_file = open(outfile2, "w")
	open_file.write(myflag + "\tPfam1_ID\tPfam1_ann\tPfam2_ID\tPfam2_ann\n")
	for myclust in sorted (DDIs.keys()):
		for item in sorted(DDIs[myclust].keys()):
			pfam1, pfam2 = item.split("\t")
			ann1 = "NA"
			ann2 = "NA"
			if pfam1 in pfam2go:
				ann1 = ";".join(sorted(pfam2go[pfam1].keys()))
			if pfam2 in pfam2go:
				ann2 = ";".join(sorted(pfam2go[pfam2].keys()))
			open_file.write(myclust + "\t" + pfam1 + "\t" + ann1 + "\t" + pfam2 + "\t" + ann2 + "\n")
		# foreach type
	# foreach cluster
	open_file.close()

	# get number
	number = {}
	for myid in anns.keys():
		mynum = len(anns[myid].keys())
		if not mynum in number:
			number[mynum] = {}
		number[mynum][myid] = ""
    # foreach DDI type

	outfile2 = re.sub(".tsv", ".info.tsv", outfile)
	open_out = open(outfile2, "w")
	open_out.write("Pfam1_ID\tPfam1_ann\tPfam2_ID\tPfam2_ann\tnumber\n")
	for mynum in sorted(number.keys(), key=int, reverse=True):
		for myid in sorted(number[mynum].keys()):	
			open_out.write(myid + "\t" + str(mynum) + "\n")
		# foreach item
	# foreah DDIs
	open_out.close()
# assign_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start ddi_DOMINE_ann.py -i " + values.ddi_annotation + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	pfams = collect_pfam_info (config.pfam_database)
	pfam2go = collect_pfam2go_info (config.pfam2go_database)
	myflag, DDI = collect_DDI_info (values.ddi_annotation)
	sys.stderr.write("Get info ......done\n")
	
	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign annotation to DDIs ......starting\n")
	assign_annotation (myflag, DDI, pfams, pfam2go, values.output)
	sys.stderr.write("\nAssign annotation to DDIs ......done\n")

	sys.stderr.write("### Finish ddi_DOMINE_ann.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

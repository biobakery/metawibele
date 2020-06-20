#!/usr/bin/env python

"""
MetaWIBELE: denovo_TM_SP module
Summary the de-novo prediction set by combining Phobius/SignalP/TMHMM

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
	from metawibele import utilities
	#import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Summary the overlapped annotation predicted by Phobius and SignalP or TMHMM
"""

def get_args (): 
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', "--phobius-signal",
	                    help='Phobius signaling prediction',
	                    required=True)
	parser.add_argument('-b', "--phobius-transmembrane",
	                    help='Phobius transmembrane prediction',
	                    required=True)
	parser.add_argument('-c', "--signlap",
	                    help='SignalP signaling prediction',
	                    required=True)
	parser.add_argument('-d', "--tmhmm",
	                    help='TMHMM transmembrane prediction',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output overlapped annotation info',
	                    required=True)
	values=parser.parse_args()

	return values
# get_args


#==============================================================
# collect annotated info
#==============================================================
def collect_annotation_info (ann_file, cluster, types):	# summary_peptide_family_annotation.tsv
	open_file = open(ann_file, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^" + utilities.PROTEIN_ID, line) or re.search("^" + utilities.PROTEIN_FAMILY_ID, line):
			title = line
			continue
		info = line.split("\t")
		myclust = info[0]
		mytype = info[1]
		types[mytype] = ""
		if not myclust in cluster:
			cluster[myclust] = {}
		cluster[myclust][mytype] = info[-1]
	# foreach line
	open_file.close()
	return title
# function collect_annnotation_info


#==============================================================
# stat overalping info between different methods
#==============================================================
def overlap_annotation(cluster, types, title, outfile):
	category_t = ["signaling", "transmembrane"]
	signaling_t = {}
	transmembrane_t = {}

	## collect types
	for mytype in types:
		if re.search("signaling", mytype) and re.search("Phobius", mytype):
			signaling_t[mytype] = ""
		if re.search("signaling", mytype) and re.search("SignalP", mytype):
			signaling_t[mytype] = ""
		if re.search("transmembrane", mytype) and re.search("Phobius", mytype):
			transmembrane_t[mytype] = ""
		if re.search("transmembrane", mytype) and re.search("TMHMM", mytype):
			transmembrane_t[mytype] = ""
	# foreach type

	## overlap info
	phobius = {}
	other = {}
	for myclust in cluster.keys():
		for mytype in cluster[myclust].keys():
			if mytype in signaling_t:
				if re.search("Phobius", mytype):
					category = "signaling"
					if not category in phobius:
						phobius[category] = {}
					phobius[category][myclust] = cluster[myclust][mytype]
				if re.search("SignalP", mytype):
					category = "signaling"
					if not category in other:
						other[category] = {}
					other[category][myclust] = cluster[myclust][mytype]
			if mytype in transmembrane_t:
				if re.search("Phobius", mytype):
					category = "transmembrane"
					if not category in phobius:
						phobius[category] = {}
					phobius[category][myclust] = cluster[myclust][mytype]
				if re.search("TMHMM", mytype):
					category = "transmembrane"
					if not category in other:
						other[category] = {}
					other[category][myclust] = cluster[myclust][mytype]
		# foreach type
	# foreach cluster
	
	# output info
	outs = {}
	for mycat in category_t:
		if mycat in phobius and mycat in other:	# extract overlaps if running both Phobius and SignalP/TMHMM
			for myclust in phobius[mycat].keys():
				myvalue1 = phobius[mycat][myclust]
				if myclust in other[mycat]:
					myvalue2 = other[mycat][myclust]
					value = myvalue1 + ";" + myvalue2
					if not mycat in outs:
						outs[mycat] = {}
					outs[mycat][myclust] = value
			# foreach cluster
		if mycat in phobius and not mycat in other: # only Phobius
			for myclust in phobius[mycat].keys():
				myvalue1 = phobius[mycat][myclust]
				if not mycat in outs:
					outs[mycat] = {}
				outs[mycat][myclust] = myvalue1
		if not mycat in phobius and mycat in other: # only SignalP/TMHMM	
			for myclust in other[mycat].keys():
				myvalue2 = other[mycat][myclust]
				if not mycat in outs:
					outs[mycat] = {}
				outs[mycat][myclust] = myvalue2
	# foreach category

	outfile1 = re.sub(".detail.tsv", ".signaling.detail.tsv", outfile)
	outfile2 = re.sub(".detail.tsv", ".transmembrane.detail.tsv", outfile)
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile2, "w")
	open_out1.write(title + "\n")
	open_out2.write(title + "\n")
	for mycat in sorted(outs.keys()):
		for myclust in sorted(outs[mycat].keys()):
			if re.search("signaling", mycat):
				mytype = "Denovo_signaling\tsignaling"
				open_out1.write(myclust + "\t" + mytype + "\t" + str(outs[mycat][myclust]) + "\n")
			if re.search("transmembrane", mycat):
				mytype = "Denovo_transmembrane\ttransmembrane"
				open_out2.write(myclust + "\t" + mytype + "\t" + str(outs[mycat][myclust]) + "\n")
		# foreach family
	# foreach category
	open_out1.close()
	open_out2.close()
# overlap_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start denovo_TM_SP.py -o " + values.output + " ####\n")

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	cluster = {}
	types = {}
	title = collect_annotation_info (values.phobius_signal, cluster, types)
	title = collect_annotation_info (values.phobius_transmembrane, cluster, types)
	title = collect_annotation_info (values.signlap, cluster, types)
	title = collect_annotation_info (values.tmhmm, cluster, types)
	sys.stderr.write("Get cluster info ......done\n")
	
	### annotation overlap ###
	sys.stderr.write("\nOverlapped annotation to peptide families ......starting\n")
	overlap_annotation(cluster, types, title, values.output)
	sys.stderr.write("\nOverlapped annotation to peptide families ......done\n")

	sys.stderr.write("### Finish denovo_TM_SP.py ####\n\n")

# end: main

if __name__ == '__main__':
	main()

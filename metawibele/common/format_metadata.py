#!/usr/bin/env python

"""
MetaWIBELE: format_metadata module
Format metadata table

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
Format metadata information
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-s', help='specify the name of study', required=True)
	parser.add_argument('-a', help='input raw metadata info', required=True)
	parser.add_argument('-b', help='input sequence info', default="none")
	parser.add_argument('-c', help='input vocabulary file to format metadata', required=True)
	parser.add_argument('-o', help='output the formated metadata file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect metadata info for samples
#==============================================================
def collect_metadata_info (study, raw_file):
	name_mapping = {"ID":"SID", "Patient":"subject", "Location":"baseline_montreal_location",
	                "Age":"consent_age", "age":"consent_age", "Location":"site_name", "HBI":"hbi", "SCCAI":"sccai", "Fecal.Calprotectin":"fecalcal"}
	meta_raw = {}
	titles = {}
	open_file = open(raw_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[info.index(item)] = item
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if study == "PRISM":
			if re.search("^[\d]+", myid):
				myid = study + "_" + myid
		myindex = 0
		while myindex < len(info):
			myitem = info[myindex]
			if myindex == 0 and study == "PRISM":
				if re.search("^[\d]+", myitem):
					myitem = study + "_" + myitem
			myname = titles[myindex]
			if myname in name_mapping:
				myname = name_mapping[myname]
			if not myid in meta_raw:
				meta_raw[myid] = {}
			meta_raw[myid][myname] = myitem
			myindex = myindex + 1
		if not "study" in meta_raw[myid]:
			meta_raw[myid]["study"] = study
	# foreach line
	open_file.close()

	return meta_raw
# collect_metadata_info


#==============================================================
# collect sequence info for samples
#==============================================================
def collect_sequence_info (study, seq_file):
	name_mapping = {"SampleName":"SID", "download_path":"data_path", "LibrarySource":"data_modality", "Run":"file_name"}
	seq = {}
	titles = {}
	title_name = {}
	open_file = open(seq_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		title_name[item] = info.index(item)
		titles[info.index(item)] = item
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if "SampleName" in title_name:
			myid = info[title_name["SampleName"]]
		if "SID" in title_name:
			myid = info[title_name["SID"]]
		myindex = 0
		while myindex < len(info):
			myitem = info[myindex]
			myname = titles[myindex]
			if myname in name_mapping:
				myname = name_mapping[myname]
			if not myid in seq:
				seq[myid] = {}
			seq[myid][myname] = myitem
			myindex = myindex + 1
		if not "study" in seq[myid]:
			seq[myid]["study"] = study
	# foreach line
	open_file.close()

	return seq
# collect_sequence_info


#==============================================================
# format info and output
#==============================================================
def format_info (meta_raw, seq, vocal_file, outfile):
	vocabulary = {}
	order = []
	titles = {}
	open_file = open(vocal_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split(",")
	for item in info:
		titles[info.index(item)] = item
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		info = line.split(",")
		myid = info[0]
		myindex = 0
		order.append(myid)
		while myindex < len(info):
			myitem = info[myindex]
			myname = titles[myindex]
			if not myid in vocabulary:
				vocabulary[myid] = {}
			vocabulary[myid][myname] = myitem
			myindex = myindex + 1
	# foreach line
	open_file.close()

	# format info
	outs = {}
	for mys in meta_raw.keys():
		skip = 0
		for mytype in vocabulary.keys():
			myvalue = "NA"
			if mytype in meta_raw[mys]:
				myvalue = meta_raw[mys][mytype]
			if mys in seq:
				if mytype in seq[mys]:
					myvalue = seq[mys][mytype]
			value_class = vocabulary[mytype]["var.class"]
			require = vocabulary[mytype]["requiredness"]
			values = vocabulary[mytype]["allowedvalues"]
			if require == "required":
				if myvalue == "NA":
					# debug
					print("No required item!\t" + mytype + "\t" + myvalue)
					skip = 1
					break
			if value_class == "numeric":
				if myvalue == "NA":
					myvalue = "NaN"
				else:
					if not myvalue.replace('.','',1).isdigit():
						# debug
						print("Not numeric value!\t" + mytype + "\t" + myvalue)
						myvalue = "NaN"
			if values != "*":
				tmp = values.split("|")
				if myvalue != "NA":
					if mytype == "sex":
						if myvalue == "f":
							myvalue = "Female"
						if myvalue == "m":
							myvalue = "Male"
					if mytype == "race":
						if myvalue == "white":
							myvalue = "White"
						if myvalue == "african_american":
							myvalue = "Black or African American"
						if myvalue == "native_american":
							myvalue = "American Indian or Alaska Native"
						if myvalue == "more_than_one":
							myvalue = "More than one race"
					if mytype == "diagnosis":
						if myvalue == "HC":
							myvalue = "nonIBD"
						if myvalue == "control":
							myvalue = "noneIBD"
					if mytype == "adult":
						if myvalue == "t" or myvalue == "true":
							myvalue = "TRUE"
						if myvalue == "f" or myvalue == "false":
							myvalue = "FALSE"
					if mytype == "pilot":
						if myvalue == "t" or myvalue == "true":
							myvalue = "TRUE"
						if myvalue == "f" or myvalue == "false":
							myvalue = "FALSE"
					if mytype == "antibiotic":
						if myvalue == "y" or myvalue == "yes":
							myvalue = "Yes"
						if myvalue == "n" or myvalue == "no":
							myvalue = "No"
					if mytype == "mesalamine":
						if myvalue == "y" or myvalue == "yes":
							myvalue = "Yes"
						if myvalue == "n" or myvalue == "no":
							myvalue = "No"
					if mytype == "steroids":
						if myvalue == "y" or myvalue == "yes":
							myvalue = "Yes"
						if myvalue == "n" or myvalue == "no":
							myvalue = "No"
					if mytype == "immunosuppressant":
						if myvalue == "y" or myvalue == "yes":
							myvalue = "Yes"
						if myvalue == "n" or myvalue == "no":
							myvalue = "No"
					if mytype == "supplement":
						if myvalue == "y" or myvalue == "yes":
							myvalue = "Yes"
						if myvalue == "n" or myvalue == "no":
							myvalue = "No"
					if mytype == "anti-diarrhoea":
						if myvalue == "y" or myvalue == "yes":
							myvalue = "Yes"
						if myvalue == "n" or myvalue == "no":
							myvalue = "No"
					if mytype == "active":
						if myvalue == "y" or myvalue == "yes":
							myvalue = "Yes"
						if myvalue == "n" or myvalue == "no":
							myvalue = "No"
					if not myvalue in tmp:
						# debug
						print("Not required value!\t" + mytype + "\t" + myvalue)
						myvalue = "NA"
			if skip == 0:
				if not mys in outs:
					outs[mys] = {}
				outs[mys][mytype] = myvalue
		# foreach type
	# foreach sample

	# output
	title = "\t".join(order)
	open_out = open(outfile, "w")
	open_out.write(title + "\n")
	for mys in sorted(outs.keys()):
		mystr = ""
		for mytype in order:
			if mytype in outs[mys]:
				myvalue = outs[mys][mytype]
			else:
				myvalue = "NA"
			if mystr == "":
				mystr = myvalue
			else:
				mystr = mystr + "\t" + myvalue
		open_out.write(mystr + "\n")
# format_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start format_metadata.py -s " + values.s + " ####\n")
	

	### format info ###
	sys.stderr.write("Format info ......starting\n")
	meta_raw = collect_metadata_info (values.s, values.a)
	if values.b == "none":
		seq = {}
	else:
		seq = collect_sequence_info (values.s, values.b)
	format_info (meta_raw, seq, values.c, values.o)
	sys.stderr.write("Format info ......done\n")
	

	sys.stderr.write("### Finish format_metadata.py ####\n\n\n")

# end: main

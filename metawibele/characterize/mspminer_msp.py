#!/usr/bin/env python

"""
MetaWIBELE: mspminer_msp module
Summary MSPs information from MSPminer

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
Summary MSPs information from MSPminer
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--annotation",
	                    help='input the uniref annotation file',
	                    required=True)
	parser.add_argument('-p', "--msp-dir",
	                    help='the output directory of MSPminer',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output summary file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect uniref annotation information
#==============================================================
def collect_uniref_info (uniref_file):
	uniref = {}
	titles = {}
	open_file = open(uniref_file, "r")
	line = open_file.readline()
	line = re.sub("\n$", "", line)
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	# foreach item
	for line in open_file:
		line = re.sub("\n$", "", line)
		if not len(line):
			continue
		info = line.split("\t")
		mymap = info[titles["map_type"]]
		myid = info[0]
		uniref[myid] = mymap
	# foreach line
	open_file.close()
	return uniref
# collect_uniref_info


#==============================================================
# collect MSPs information
#==============================================================
def collect_msp_info (msp_dir, uniref, outfile):
	msp = {}
	msp_stat = {}
	for mymsp in os.listdir(msp_dir):
		mymsp_dir = os.path.join(msp_dir, mymsp)
		if not os.path.isdir(mymsp_dir):
			continue
		myfile = mymsp_dir + "/modules.tsv"
		open_file = open(myfile, "r")
		for line in open_file:
			line = line.strip()
			if not len(line):
				continue
			if re.search("^module_name", line):
				continue
			info = line.split("\t")
			module_name = info[0]
			gene_id = info[1]
			gene_name = info[2]
			mytype = module_name
			mytype = re.sub("_module[\S]+", "", mytype)
			if not mymsp in msp:
				msp[mymsp] = {}
			if not module_name in msp[mymsp]:
				msp[mymsp][module_name] = {}
			msp[mymsp][module_name][gene_name] = mytype + "\t" + module_name + "\t" + gene_id + "\t" + gene_name
			if not mymsp in msp_stat:
				msp_stat[mymsp] = {}
			if not mytype in msp_stat[mymsp]:
				msp_stat[mymsp][mytype] = {}
			msp_stat[mymsp][mytype][gene_name] = ""
		# foreach line
		open_file.close()
	# foreach MSP

	# output info
	open_out = open(outfile, "w")
	open_out.write("msp_name\tclass\tmodule_name\tgene_id\tgene_name\tcmp_type\n")
	for mymsp in sorted(msp.keys()):
		for mymodule in sorted(msp[mymsp].keys()):
			for mygene in sorted(msp[mymsp][mymodule].keys()):
				mymap = "NA"
				if mygene in uniref:
					mymap = uniref[mygene]
				open_out.write(mymsp + "\t" + msp[mymsp][mymodule][mygene] + "\t" + mymap + "\n")
	# foreach MSP
	open_out.close()

	outfile1 = re.sub(".tsv", ".stat.tsv", outfile)
	open_out = open(outfile1, "w")
	open_out.write("msp_name\t# genes\t# core genes\t# accessory genes\t# shared core genes\t# shared accessory genes\n")
	for mymsp in sorted(msp_stat.keys()):
		mytotal = 0
		mycore = 0
		myacc = 0
		myshared_core = 0
		myshared_acc = 0
		if "core" in msp_stat[mymsp]:
			mycore = len(msp_stat[mymsp]["core"])
		if "accessory" in msp_stat[mymsp]:
			myacc = len(msp_stat[mymsp]["accessory"])
		if "shared_core" in msp_stat[mymsp]:
			myshared_core = len(msp_stat[mymsp]["shared_core"])
		if "shared_accessory" in msp_stat[mymsp]:
			myshared_acc = len(msp_stat[mymsp]["shared_accessory"])
		mytotal = mycore + myacc + myshared_core + myshared_acc
		open_out.write(mymsp + "\t" + str(mytotal) + "\t" + str(mycore) + "\t" + str(myacc) + "\t" + str(myshared_core) + "\t" + str(myshared_acc) + "\n")
	# foreach MSP
	open_out.close()
# func: collect_msp_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start mspminer_msp.py -a " + values.annotation + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	uniref = collect_uniref_info (values.annotation)
	collect_msp_info (values.msp_dir, uniref, values.output)
	sys.stderr.write("Get info ......done\n")


	sys.stderr.write("### Finish mspminer_msp.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

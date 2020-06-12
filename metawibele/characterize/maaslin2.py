#!/usr/bin/env python

"""
MetaWIBELE: maaslin2 module
Run maaslin2 to associate features with metadata

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
Run MaAsLin2
"""
	
def get_args():	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', "--feature-table",
	                    help='input feature table, features as columns and samples as rows',
	                    required=True)
	parser.add_argument('-m', "--metadata-table",
	                    help='input metadata file',
	                    required=True)
	parser.add_argument('-n', "--split-num",
	                    help='number of split files',
	                    default=1)
	parser.add_argument('-w', "--workdir",
	                    help='workding directory',
	                    default="none")
	parser.add_argument('-o', "--output",
	                    help='output file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect info
#==============================================================
def collect_metadata (meta_file, spe_type, meta, outfile):
	titles = {}
	open_file = open(meta_file, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	open_out.write(line)
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = re.sub("\n$", "", line)
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if spe_type in titles:
			myindex = titles[spe_type]
			mymeta = info[myindex]
			if mymeta == meta:
				open_out.write(line + "\n")
	# foreach line
	open_file.close()
	open_out.close()
# collect_metadata


#==============================================================
# run MAasLin2
#==============================================================
def run_maaslin2 (feature, metadata, split_num, workdir, output, outfile):
	myexe = config.maaslin2_cmmd
	myopt = " ".join([str(i) for i in (config.maaslin2_cmmd_opts)])
	mytranspos = config.transpose_cmmd
	mypcl = config.pcl_utils
	myutils = config.maaslin2_utils

	if split_num == 1:
		mycmd = "Rscript " + myexe + " " + feature + " " + metadata + " " + output + " " + myopt
		print(mycmd)
		os.system(mycmd)
	else:
		mytrans = re.sub(".tsv", ".pcl", feature)
		mym = re.search("([^\/]+)$", mytrans)
		mybase = mym.group(1)
		mybase = re.sub(".pcl", ".split", mybase)
		if not os.path.isfile(mytrans):
			os.system(mytranspos + " < " + feature + " > " + mytrans)
		# split files
		mycmd = "Rscript " + myutils + " " + "split" + " " + mytrans + " " + mybase + " " + mypcl + " " + str(split_num)
		print(mycmd)
		os.system(mycmd)
		mylist = re.sub(".tsv", ".pcl.list", feature)
		os.system("ls " + mybase + "*.pcl > " + mylist)
		os.system("mkdir -p " + output)
		os.system("mv " + mybase + "*.pcl" + " " + output)
		myout = workdir + "/" + outfile
		os.system("rm -f " + myout)
		open_file = open(mylist, "r")
		for i in open_file:
			i = i.strip()
			if not len(i):
				continue
			myfeature = re.sub(".pcl", ".tsv", i)
			myinput = output + "/" + myfeature
			mycmd = mytranspos + " < " + output + "/" + i + " > " + myinput
			print(mycmd)
			os.system(mycmd)
			myoutput = re.sub(".pcl", "", i)
			myoutput = output + "/" + myoutput
			mycmd = "Rscript " + myexe + " " + myinput + " " + metadata + " " + myoutput + " " + myopt
			print(mycmd)
			os.system(mycmd)
			if not os.path.isfile(myout):
				os.system("less " + myoutput + "/all_results.tsv > " + myout)
			else:
				os.system("sed -e 1d " + myoutput + "/all_results.tsv >> " + myout)
			
			# delete intermediate split files
			mypcl = output + "/" + i
			os.system("rm -f " + mypcl)
			os.system("rm -f " + myinput)
		# foreach split file
		open_file.close()
	# if split files

# run_maaslin2


def fdr_correction (output, outfile):
	myin = output + "/" + outfile
	myout = re.sub(".tsv", ".fdr_correction.tsv", myin)
	mypcl = config.pcl_utils
	myutils = config.maaslin2_utils
	mycmd = "Rscript " + myutils + " " + "correct" + " " + myin + " " + myout + " " + mypcl + " 0 "
	print(mycmd)
	os.system(mycmd)
	myout1 = myout
	myout2 = re.sub(".tsv", ".correct_per_variable.tsv", myout)
	myout3 = re.sub(".tsv", ".correct_per_level.tsv", myout)
	return myout1, myout2, myout3
# fdr_correction


#==============================================================
###########  Main processing ############
#==============================================================
def main():

	### get arguments ###

	values = get_args()
	if values.workdir == "none":
		values.workdir = config.maaslin2_dir

	try: 
		values.split_num = int(values.split_num)
	except ValueError:
		sys.exit("Please specify valid number for spliting file")

	sys.stderr.write("### Start maaslin2.py -i " + values.feature_table + " ####\n")
	
	### Run MaAsLin2 info ###
	if config.nested_effects != "none":
		out1 = []
		out2 = []
		out3 = []
		tmp1 = config.nested_effects.split(":")
		mytype = tmp1[0]
		tmp2 = tmp1[1].split(",")
		for item in tmp2:
			#output = config.maaslin2_dir + "/" + item
			output = values.workdir + "/" + item
			mymeta_file = re.sub(".tsv", "." + item + ".tsv", values.metadata_table)
			collect_metadata (values.metadata_table, mytype, item, mymeta_file)
			outfile = re.sub(".tsv", "." + item + ".tsv", values.output)
			sys.stderr.write("Run MaAsLin2 ......" + mymeta_file + "\n")
			run_maaslin2 (values.feature_table, mymeta_file, values.split_num, values.workdir, output, outfile)
			sys.stderr.write("Run MaAsLin2 ......done\n")
			sys.stderr.write("FDR correction ......" + outfile + "\n")
			myout1, myout2, myout3 = fdr_correction (values.workdir, outfile)
			sys.stderr.write("FDR correction ......done\n")
			out1.append(myout1)
			out2.append(myout2)
			out3.append(myout3)
		# foreach nested effect

		# collect all results
		#outfile = config.maaslin2_dir + "/" + values.output
		outfile = values.workdir + "/" + values.output
		outfile1 = re.sub(".tsv", ".fdr_correction.tsv", outfile)
		outfile2 = re.sub(".tsv", ".fdr_correction.correct_per_variable.tsv", outfile)
		outfile3 = re.sub(".tsv", ".fdr_correction.correct_per_level.tsv", outfile)
		os.system("less " + out1[0] + " > " + outfile1)
		os.system("less " + out2[0] + " > " + outfile2)
		os.system("less " + out3[0] + " > " + outfile3)
		myindex = 1
		while myindex < len(out1):
			os.system("sed -e 1d " + out1[myindex] + " >> " + outfile1)	
			os.system("sed -e 1d " + out2[myindex] + " >> " + outfile2)	
			os.system("sed -e 1d " + out3[myindex] + " >> " + outfile3)	
			myindex = myindex + 1
	# if nested effect
	else:
		sys.stderr.write("Run MaAsLin2 ......starting\n")
		run_maaslin2 (values.feature_table, values.metadata_table, values.split_num, values.workdir, values.workdir, values.output)
		sys.stderr.write("Run MaAsLin2 ......done\n")

		sys.stderr.write("FDR correction ......starting\n")
		fdr_correction (values.workdir, values.output)
		sys.stderr.write("FDR correction ......done\n")
	
	sys.stderr.write("### Finish maaslin2.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

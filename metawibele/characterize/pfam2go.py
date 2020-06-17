#!/usr/bin/env python

"""
MetaWIBELE: pfam2go module
Summary GO annotation info based on Pfam2GO file

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
Summary GO annotation info based on Pfam2GO file
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--pfam-annotation",
	                    help='input Pfam annotation file',
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
def collect_pfam2go_info (annfile):	# Pfam2GO.txt 
	pfam2go = {}
	for line in utils.gzip_bzip2_biom_open_readlines (annfile):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		info = line.split("\t")
		mypfam = info[0]
		mygo = info[1]
		mydec = info[2]
		mytype = info[3]
		if mytype == "process":
			mytype = "BP"
		if mytype == "function":
			mytype = "MF"
		if mytype == "component":
			mytype = "CC"
		if not mypfam in pfam2go:
			pfam2go[mypfam] = {}
		pfam2go[mypfam][mygo] = mydec + "\t" + mytype
	# foreach line
	
	return pfam2go
# function collect_pfam2go_info


def collect_pfam_info (pfamfile):	# summary_Pfam_peptide.detail.tsv
	pfams = {}
	open_file = open(pfamfile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myflag = info[0]
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		mym = re.search("([^\_]+)$", info[1])
		mytype = mym.group(1)
		if mytype != "PfamDomain":
			continue
		pfams[myid] = info[2] + "\t" + info[-1]
	# foreach line
	open_file.close()
	return myflag, pfams
#collect_pfam_info


#==============================================================
# assign annotation
#==============================================================
def assign_annotation (id_flag, pfams, pfam2go, outfile):
	#outfile1 = re.sub(".tsv", ".ontologizer.association.tsv", outfile)
	#outfile2 = re.sub(".tsv", ".genelist.tsv", outfile)
	outfile1 = re.sub(".tsv", ".split.tsv", outfile)
	open_file = open(outfile, "w")
	open_file1 = open(outfile1, "w")
	open_file.write(id_flag + "\ttype\tdetail\tdescription\tconsistency\n")
	open_file1.write(id_flag + "\ttype\tdetail\tdescription\tconsistency\n")
	for myclust in sorted (pfams.keys()):
		mypfam, mycons = pfams[myclust].split("\t")
		tmp = mypfam.split(";")
		tmp_cons = mycons.split(";")
		mystr1_1 = {}
		mystr2_1 = {}
		mystr3_1 = {}
		myindex = 0
		while myindex < len(tmp):
			item = tmp[myindex]
			mycon = tmp_cons[myindex]
			if item in pfam2go:
				for mygo in sorted(pfam2go[item].keys()):
					tmp2 = pfam2go[item][mygo].split("\t")
					mytype = "Pfam2GO_GO(" + tmp2[-1] + ")"
					mystr = myclust + "\t" + mytype + "\t" + mygo + "\t" + tmp2[0] + "\t" + mycon
					open_file1.write(mystr + "\n")
					if tmp2[-1] == "BP":
						mystr1_1[mygo] = tmp2[0] + "\t" + mycon
					if tmp2[-1] == "CC":
						mystr2_1[mygo] = tmp2[0] + "\t" + mycon
					if tmp2[-1] == "MF":
						mystr3_1[mygo] = tmp2[0] + "\t" + mycon
			# if go info
			myindex = myindex + 1
		# foreach pfam

		if len(mystr1_1.keys()) > 0:
			mygo = ";".join(sorted(mystr1_1.keys()))
			myann = "NA"
			mycon = "NA"
			for i in sorted(mystr1_1.keys()):
				tmp = mystr1_1[i].split("\t")
				if myann == "NA":
					myann = tmp[0]
					mycon = tmp[1]
				else:
					myann = myann + ";" + tmp[0]
					mycon = mycon + ";" + tmp[1]
			open_file.write(myclust + "\t" + "Pfam2GO_GO(BP)" + "\t" + mygo + "\t" + myann + "\t" + mycon + "\n")
		if len(mystr2_1.keys()) > 0:
			mygo = ";".join(sorted(mystr2_1.keys()))
			myann = "NA"
			mycon = "NA"
			for i in sorted(mystr2_1.keys()):
				tmp = mystr2_1[i].split("\t")
				if myann == "NA":
					myann = tmp[0]
					mycon = tmp[1]
				else:
					myann = myann + ";" + tmp[0]
					mycon = mycon + ";" + tmp[1]
			open_file.write(myclust + "\t" + "Pfam2GO_GO(CC)" + "\t" + mygo + "\t" + myann + "\t" + mycon + "\n")
		if len(mystr3_1.keys()) > 0:
			mygo = ";".join(sorted(mystr3_1.keys()))
			myann = "NA"
			mycon = "NA"
			for i in sorted(mystr3_1.keys()):
				tmp = mystr3_1[i].split("\t")
				if myann == "NA":
					myann = tmp[0]
					mycon = tmp[1]
				else:
					myann = myann + ";" + tmp[0]
					mycon = mycon + ";" + tmp[1]
			open_file.write(myclust + "\t" + "Pfam2GO_GO(MF)" + "\t" + mygo + "\t" + myann + "\t" + mycon + "\n")
	# foreach cluster
	open_file.close()
	open_file1.close()
# assign_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start pfam2go.py -i " + values.pfam_annotation + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	id_flag, pfams = collect_pfam_info (values.pfam_annotation)
	pfam2go = collect_pfam2go_info (config.pfam2go_database)
	sys.stderr.write("Get info ......done\n")
	
	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign annotation ......starting\n")
	assign_annotation (id_flag, pfams, pfam2go, values.output)
	sys.stderr.write("\nAssign annotation ......done\n")

	sys.stderr.write("### Finish pfam2go.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

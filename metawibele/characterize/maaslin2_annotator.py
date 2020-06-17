#!/usr/bin/env python

"""
MetaWIBELE: maaslin2_annotator module
Summary maaslin2 annotation

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
import math

try:
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Summary association annotation for protein families
"""

def get_args():

	parser=argparse.ArgumentParser()
	parser.add_argument('-s', "--stat",
	                    help='input statistical summary file',
	                    required=True)
	parser.add_argument('-t', "--type",
	                    help='specify the type of abundance table, e.g. MaAsLin2_DA',
	                    required=True)
	parser.add_argument('-e', "--effect-size",
	                    help='specify the item name indicating effect size', 
						required=True)
	parser.add_argument('-o', "--output",
	                    help='output annotated file',
	                    required=True)
	values = parser.parse_args()

	return values
# get_args


#==============================================================
# get annotation about differential abundance info
#==============================================================
def stat_annotation (statfile, mytype, effect_name, outfile):
	titles = {}
	stat = []
	open_file = open(statfile, "r")
	open_out = open(outfile, "w")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_FAMILY_ID, line):
			for item in info:
				titles[item] = info.index(item)
			mytitle = "\t".join(info[1:])
			open_out.write(utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\t" + mytitle + "\n")
			continue
		myid = info[titles[utilities.PROTEIN_FAMILY_ID]]
		myeffect = info[titles["coef"]]
		if effect_name in titles:
			myeffect = info[titles[effect_name]]
		mycmp = "NA"
		if "cmp_type" in titles:
			mycmp = info[titles["cmp_type"]]
		mydetail = "NA"
		if myeffect != "NA" and myeffect != "NaN" and myeffect != "nan":
			if float(myeffect) < 0:
				mydetail = "down"
			else:
				if float(myeffect) > 0:
					mydetail = "up"
				else:
					mydetail = "NA"
		if mycmp != "NA":
			mydetail = mydetail + "(" + mycmp + ")"
		mystat = "\t".join(info[1:])
		open_out.write (myid + "\t" + mytype + "\t" + mydetail + "\t" + mystat + "\n")
	# foreach line
	open_file.close()
	open_out.close()

# func: stat_annotation



#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args()

	sys.stderr.write("### Start maaslin2_annotator.py -s " + values.stat + " ####\n")

	### get info
	sys.stderr.write("\nGet DA annotation info ......starting\n")
	values.effect_size = re.sub("mean_log", "mean(log)", values.effect_size)
	values.effect_size = re.sub("log_FC", "log(FC)", values.effect_size)
	stat_annotation (values.stat, values.type, values.effect_size, values.output)
	sys.stderr.write("Get DA annotation info ......done\n")

	sys.stderr.write("### Finish maaslin2_annotator.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

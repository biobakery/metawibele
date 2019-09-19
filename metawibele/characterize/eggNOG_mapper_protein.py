#!/usr/bin/env python

"""
MetaWIBELE: eggNOG_mapper_protein module
Extract predicted proteins based on eggNOG-Mapper results

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
Extract predicted peptides based on eggNOG-Mapper results
"""

def get_args ():
	parser = argparse.ArgumentParser()
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
						default="emapper.annotations")
	parser.add_argument('-p', "--path",
	                    help='input the path of eggNOG-Mapper file',
	                    required=True)
	values = parser.parse_args()	
	return values
# get_args


#==============================================================
# collect eggNOG-Mapper info
#==============================================================
def extract_emapper_info (GOs, extension, emapper_path):
	filelist = utilities.find_files(emapper_path, extension, None)
	for myfile in filelist:
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			print ("OK!\t" + myfile)
			myout1 = re.sub(".emapper.annotations", ".emapper.GO.tsv", myfile)
			myout2 = re.sub(".emapper.annotations", ".emapper.KEGG-KOs.tsv", myfile)
			myout3 = re.sub(".emapper.annotations", ".emapper.BiGG.tsv", myfile)
			myout4 = re.sub(".emapper.annotations", ".emapper.COG.tsv", myfile)
			open_out1 = open(myout1, "w")
			open_out2 = open(myout2, "w")
			open_out3 = open(myout3, "w")
			open_out4 = open(myout4, "w")

			open_file = open(myfile, "r")
			open_out1.write("seqID\tGO_terms\tannotation\ttype\n")
			open_out2.write("seqID\tKEGG-KOs\tannotation\ttype\n")
			open_out3.write("seqID\tBiGG_reactions\tannotation\ttype\n")
			open_out4.write("seqID\tCOG\tannotation\ttype\n")
			titles = {}
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^#", line):
					if re.search("^#query_name", line):
						info = line.split("\t")
						for item in info:
							titles[item] = info.index(item)
					continue
				info = line.split("\t")
				myid = info[titles["#query_name"]]
				go = info[titles["GO_terms"]]
				if go == "":
					go = "NA"
				kegg = info[titles["KEGG_KOs"]]
				if kegg == "":
					kegg = "NA"
				bigg = info[titles["BiGG_reactions"]]
				if bigg == "":
					bigg = "NA"
				ogs = info[titles["OGs"]]
				if ogs == "":
					ogs = "NA"
				best_og = info[titles["bestOG|evalue|score"]]
				tmp = best_og.split("|")
				if tmp == 0:
					best_og = "NA"
				else:
					best_og = tmp[0]
				og_cat = info[titles["COG cat"]]
				if og_cat == "":
					og_cat = "NA"
				og_ann = info[titles["eggNOG annot"]]
				if og_ann == "":
					og_ann = "NA"
				
				# GO
				tmp = go.split(",")
				for item in tmp:
					if item in GOs:
						myacc, myann = GOs[item].split("\t")
					else:
						myacc = item
						myann = "GO"
					if myacc != "NA[NA]" and myacc != "NA":
						mydesc = "NA"
						if re.search("\[GO", myacc):
							mym = re.search("([^\[]+)\[([^\]]+)\]", myacc)
							mydesc = mym.group(1)
							myacc = mym.group(2)
						open_out1.write(myid + "\t" + myacc + "\t" + mydesc + "\t" + myann + "\n")
				
				# KEGG
				tmp = kegg.split(",")
				for item in tmp:
					myann = "KEGG-KOs"
					if item != "NA":
						open_out2.write(myid + "\t" + item + "\tNA\t" + myann + "\n")
				
				# BiGG
				tmp = bigg.split(",")
				for item in tmp:
					myann = "BiGG"
					if item != "NA":
						open_out3.write(myid + "\t" + item + "\tNA\t" + myann + "\n")
				
				# COG
				if best_og != "NA":
					open_out4.write(myid + "\t" + best_og + "\t" + og_ann + "#" + og_cat + "#" + ogs + "\tCOG" + "\n")
			# foreach line	
			open_out1.close()
			open_out2.close()
			open_out3.close()
			open_out4.close()
			open_file.close()
		# if file exist
	# foreach samplelist
# function extract_interproscan_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start eggNOG_mapper_protein.py -e " + values.extension + " ####\n")
	
	### collect info ###
	sys.stderr.write("Get eggNOG-Mapper info ......starting\n")
	GOs = utilities.collect_GO_info (config.go_database)
	extract_emapper_info (GOs, values.extension, values.path)
	sys.stderr.write("Get eggNOG-Mapper info ......done\n")

	sys.stderr.write("### Finish eggNOG_mapper_protein.py ####\n\n\n")

# end: main

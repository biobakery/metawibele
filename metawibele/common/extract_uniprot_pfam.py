#!/usr/bin/env python
##########################################################################
# Function: Extract the Pfam and UniProt mapping info (especially for Human domain) 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 09/15/2018 
##########################################################################
import sys
import os
import re
import argparse


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Extract the Pfam and UniProt mapping info
"""

def get_args ():
	parser=argparse.ArgumentParser()
	parser.add_argument('-a', help='the input UniProt annotation file', required=True)
	parser.add_argument('-o', help='the output file for Pfam ID and corresponding UniProt info', required=True)
	values=parser.parse_args()
	return values
# get_args


#==============================================================
# Extract annotation info
#==============================================================
def extract_annotation_info (seqfile):	# uniprot_annotation
	titles = {}
	ann1 = {}
	ann2 = {}
	ann3 = {}
	ann4 = {}
	ann5 = {}
	open_file = open (seqfile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 0
	items = ["Entry name", "Gene names", "Protein names", "Organism", "NCBI_TaxID", "Cross-reference (Pfam)"]
	while myindex < len(info):
		item = info[myindex]
		if item == "Taxonomic lineage IDs":
			item = "NCBI_TaxID"
		if item in items:
			titles[myindex] = item
		myindex = myindex + 1
	# foreach item
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myindex = 1
		myid = "NA"
		org = "NA"
		taxa = "NA"
		gene = "NA"
		protein = "NA"
		pfam_item = "NA"
		while myindex < len(info):
			item = info[myindex]
			if myindex in titles:
				# UniProtKB-ID
				if titles[myindex] == "Entry name" and item != "NA" and item != "":
					myid = item
				# Gene name
				if titles[myindex] == "Gene names" and item != "NA" and item != "":
					gene = item
					gene = re.sub(";$", "", gene)
					gene = re.sub("\{[\S\s]+\}", "", gene)
					gene = re.sub("\s+", ";", gene)
				# Protein name
				if titles[myindex] == "Protein names" and item != "NA" and item != "":
					protein = item
					#protein = re.sub("\{ECO[\S\s]+\}", "", protein)
				# Organism
				if titles[myindex] == "Organism" and item != "NA" and item != "":
					org = item
				# taxanomy
				if titles[myindex] == "NCBI_TaxID" and item != "NA" and item != "":
					taxa = item
				# Pfam
				if titles[myindex] == "Cross-reference (Pfam)" and item != "NA" and item != "":
					pfam_item = item
			# if items
			myindex = myindex + 1
		# foreach item
		if pfam_item != "NA":
			pfam_item = re.sub("\s+", "", pfam_item)
			items = pfam_item.split(";")
			for pfam in items:
				if pfam == "":
					continue
				if not pfam in ann1:
					ann1[pfam] = {}
				ann1[pfam][org] = ""
				if not pfam in ann2:
					ann2[pfam] = {}
				ann2[pfam][myid] = ""
				if not pfam in ann3:
					ann3[pfam] = {}
				ann3[pfam][protein] = ""
				if not pfam in ann4:
					ann4[pfam] = {}
				ann4[pfam][gene] = ""
				if not pfam in ann5:
					ann5[pfam] = {}
				ann5[pfam][taxa] = ""
			# foreach pfam
		# foreach item
	# foreach line
	open_file.close()
	return ann1, ann2, ann3, ann4, ann5
# extract_annotation_info


#==============================================================
# output info
#==============================================================
def output_info (ann1, ann2, ann3, ann4, ann5, outfile):
	open_file = open(outfile, "w")
	open_file.write("Pfam\tOrganism\tNCBI_TaxID\tUniProtKB\tProtein_name\tGene_name\n")
	for mypfam in sorted(ann1.keys()):
		org = ";".join(sorted(ann1[mypfam].keys()))
		taxa = ";".join(sorted(ann5[mypfam].keys()))
		uniprot = ";".join(sorted(ann2[mypfam].keys()))
		protein = ";".join(sorted(ann3[mypfam].keys()))
		gene = ";".join(sorted(ann4[mypfam].keys()))
		org = re.sub("NA;", "", org)
		org = re.sub(";NA", "", org)
		taxa = re.sub("NA;", "", taxa)
		taxa = re.sub(";NA", "", taxa)
		uniprot = re.sub("NA;", "", uniprot)
		uniprot = re.sub(";NA", "", uniprot)
		protein = re.sub("NA;", "", protein)
		protein = re.sub(";NA", "", protein)
		gene = re.sub("NA;", "", gene)
		gene = re.sub(";NA", "", gene)
		open_file.write(mypfam + "\t" + org + "\t" + taxa + "\t" + uniprot + "\t" + protein + "\t" + gene + "\n")
	# foreach Pfam
	open_file.close()
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():

	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start extract_uniprot_pfam.py -a " + values.a + " ####\n")
	
	### Extraction ###
	sys.stderr.write("Extract info......starting\n")
	ann1, ann2, ann3, ann4, ann5 = extract_annotation_info (values.a)
	sys.stderr.write("Extract info......done\n")
	
	### Output ###
	sys.stderr.write("\nOutput info......starting\n")
	output_info (ann1, ann2, ann3, ann4, ann5, values.o)
	sys.stderr.write("\nOutput info......done\n")
	
	sys.stderr.write("### Finish  extract_uniprot_pfam.py  ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python

"""
MetaWIBELE: summary_all_annotation module
Summary functional and taxonomic annotations

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
Summary functional and taxonomic annotation
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--annotation",
	                    help='input functional annotation file',
	                    required=True)
	parser.add_argument('-t', "--taxonomy",
	                    help='input taxonomic annotation file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output summarized total annotations file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect functional annotation info
#==============================================================
def collect_function_info (ann_file):  
	note = {}
	anns = {}
	titles = {}
	dna = {}
	rna = {}
	open_file = open(ann_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	id_flag = info[0]
	for item in info:
		titles[item] = info.index(item)
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		mytype = info[titles["type"]]
		myann = info[titles["detail"]]
		mynote = info[titles["note"]]
		note[info[0]] = mynote
		if mytype == "DNA_abundance":
			if not myid in dna:
				dna[myid] = {}
			dna[myid]["DNA_abundance"] = myann
		if mytype == "DNA_prevalence":
			if not myid in dna:
				dna[myid] = {}
			dna[myid]["DNA_prevalence"] = myann
		if mytype == "RNA_abundance":
			if not myid in rna:
				rna[myid] = {}
			rna[myid]["RNA_abundance"] = myann
		if mytype == "RNA_prevalence":
			if not myid in rna:
				rna[myid] = {}
			rna[myid]["RNA_prevalence"] = myann
		if mytype == "UniRef90_unknown" or mytype == "UniRef90_uncharacterized":
			continue
		if not myid in anns:
			anns[myid] = mytype + "=" + myann
		else:
			anns[myid] = anns[myid] + "#" + mytype + "=" + myann
	# foreach line
	open_file.close()
	return anns, dna, rna, note
# function collect_function_info


#==============================================================
# combine functional and taxonomic annotation info
#==============================================================
def combine_annotation (annotation, dna, rna, note, taxa_file, outfile):
	titles = {}
	open_file = open(taxa_file, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	rna_flag = 0
	if len(rna.keys()) > 0:
		rna_flag = 1
	if rna_flag == 1:
		title = info[0] + "\t" + info[1] + "\tmap_type\tunirefID\tUniProtKB\tProtein_names\tmsp_name\tmsp_taxa_name\tmsp_taxa_id\ttaxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\tDNA_abundance\tDNA_prevalence\tRNA_abundance\tRNA_prevalence\tannotation\tnote"
	else:
		title = info[0] + "\t" + info[1] + "\tmap_type\tunirefID\tUniProtKB\tProtein_names\tmsp_name\tmsp_taxa_name\tmsp_taxa_id\ttaxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\tDNA_abundance\tDNA_prevalence\tannotation\tnote"
	open_out.write(title + "\n")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		mystudy = info[1]
		mymap = info[titles["map_type"]]
		myuniref = info[titles["unirefID"]]
		myunprot = info[titles["UniProtKB"]]
		mydesc = info[titles["detail"]]
		mymsp = info[titles["msp_name"]]
		mymsp_taxa = info[titles["msp_taxa_name"]]
		mymsp_taxa_id = info[titles["msp_taxa_id"]]
		taxa_id = info[titles["taxa_id"]]
		taxa_name = info[titles["taxa_name"]]
		taxa_rank = info[titles["taxa_rank"]]
		taxa_lineage = info[titles["taxa_lineage"]]
		mynote = info[titles["note"]]
		mydna_a = "NA"
		mydna_p = "NA"
		myrna_a = "NA"
		myrna_p = "NA"
		if myid in dna:
			if "DNA_abundance" in dna[myid]:
				mydna_a = dna[myid]["DNA_abundance"]
			if "DNA_prevalence" in dna[myid]:
				mydna_p = dna[myid]["DNA_prevalence"]
		if myid in rna:
			if "RNA_abundance" in rna[myid]:
				myrna_a = rna[myid]["RNA_abundance"]
			if "RNA_prevalence" in rna[myid]:
				myrna_p = rna[myid]["RNA_prevalence"]
		if myid in note:
			mynote1 = note[myid]
			tmp1 = mynote1.split(";")
			tmp2 = mynote.split(";")
			note_tmp = {}
			for item in tmp1:
				note_tmp[item] = ""
			for item in tmp2:
				note_tmp[item] = ""
			mynote = ";".join(sorted(note_tmp))
		myann = "NA"
		if myid in annotation:
			myann = annotation[myid]
		if rna_flag == 1:
			mystr = myid + "\t" + mystudy + "\t" + mymap + "\t" + myuniref + "\t" + myunprot + "\t" + mydesc + "\t" + mymsp + "\t" + mymsp_taxa + "\t" + mymsp_taxa_id + "\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage + "\t" + mydna_a + "\t" + mydna_p + "\t" + myrna_a + "\t" + myrna_p + "\t" + myann + "\t" + mynote
		else:
			mystr = myid + "\t" + mystudy + "\t" + mymap + "\t" + myuniref + "\t" + myunprot + "\t" + mydesc + "\t" + mymsp + "\t" + mymsp_taxa + "\t" + mymsp_taxa_id + "\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage + "\t" + mydna_a + "\t" + mydna_p + "\t" + myann + "\t" + mynote
		open_out.write(mystr + "\n")
	# foreach cluster
	open_file.close()
	open_out.close()

# combine_annotation 


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start summary_all_annotation.py -a " + values.annotation + " ####\n")
	

	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	anns, dna, rna, note = collect_function_info (values.annotation)
	sys.stderr.write("Get annotation info ......done\n")

	### combine functional and taxonomic annotation
	sys.stderr.write("\nCombine annotation......starting\n")
	combine_annotation (anns, dna, rna, note, values.taxonomy, values.output)
	sys.stderr.write("\nCombine annotation......done\n")

	sys.stderr.write("### Finish summary_all_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

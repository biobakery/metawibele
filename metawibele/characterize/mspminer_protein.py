#!/usr/bin/env python

"""
MetaWIBELE: mspminer_protein module
Extract annotation of MSPs

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
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Extract MSP annotation for each protein 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--annotation",
	                    help='input the MSP taxonomy annotation file',
	                    required=True)
	parser.add_argument('-c', "--type",
	                    help='specify the class of genes used for taxonomy annotation',
	                    choices=["all", "core"],
	                    default="all")
	parser.add_argument('-m', "--msp",
	                    help='input the MSP annotation information file',
	                    required=True)
	parser.add_argument('-g', "--gene",
	                    help='input the gene file for MSP annotation',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output summary file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect taxonomy annotation info of MSPs
#==============================================================
def collect_taxonomy_info (ann_file, class_type): 
	anns = {}
	titles = {}
	open_file = open(ann_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("msp_name\t", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		# if title
		mymsp = info[titles["msp_name"]]
		myclass = info[titles["class"]]
		num_genes = info[titles["num_genes"]]
		num_unclassified_genes = info[titles["num_unclassified_genes"]]
		num_hit_genes = info[titles["num_hit_genes"]]
		taxa_id = info[titles["taxa_id"]]
		taxa_name = info[titles["taxa_name"]]
		taxa_rank = info[titles["taxa_rank"]]
		taxa_lineage = info[titles["taxa_lineage"]]
		if myclass != class_type:
			continue
		if not mymsp in anns:
			anns[mymsp] = "num_genes(" + num_genes + ")," + "num_unclassified_genes(" + num_unclassified_genes + ")," + "num_hit_genes(" + num_hit_genes + ")\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage
	# foreach line
	open_file.close()
	return anns
# collect_taxonomy_info


#==============================================================
# collect gene info of MSPs
#==============================================================
def collect_msp_info (msp_file): 
	msp = {}
	titles = {}
	open_file = open(msp_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("msp_name\t", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		# title line
		mymsp = info[titles["msp_name"]]
		myid = info[titles["gene_name"]]
		myclass = info[titles["msp_name"]] + "__" + info[titles["class"]]
		#myclass = info[titles["class"]]
		mymodule = info[titles["msp_name"]] + "__" + info[titles["module_name"]]
		msp[myid] = mymsp + "\n" + myclass + "\t" + mymodule
	# foreach line
	open_file.close()
	return msp
# foreach collect_msp_info


#==============================================================
# assign MSP annotation info to gene
#==============================================================
def assign_msp (gene_file, msp, anns, outfile):  
	genes = {}
	titles = {}
	taxa_title = ""
	open_file = open(gene_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search(utilities.PROTEIN_ID, line):
			taxa_title = line
			for item in info:
				titles[item] = info.index(item)
			continue
		genes[info[0]] = line
	# foreach line
	open_file.close()

	open_out = open(outfile, "w")
	open_out.write(utilities.PROTEIN_ID + "\tmsp_name\tgene_class\tmodule_name\tmsp_description\ttaxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\n")
	for myid in sorted(genes.keys()):
		if myid in msp:	# grouped into MSP
			mymsp, myinfo = msp[myid].split("\n")
			myann = "NA\tNA\tNA\tUnclassified\tNA"
			if mymsp in anns:
				myann = anns[mymsp]
			else:
				# debug
				print("No MSP annotation info!\t" + mymsp)
			mystr = myid + "\t" + mymsp + "\t" + myinfo + "\t" + myann
			open_out.write(mystr + "\n")
		else:	# not grouped into MSP
			open_out.write(myid + "\t" + "msp_unknown" + "\tNA\tNA" + "\t" + "NA\tNA\tNA\tUnclassified\tNA" + "\n")
	# foreach gene
	open_out.close()
	
	outfile1 = re.sub(".tsv", ".taxonomy.tsv", outfile)
	open_out = open(outfile1, "w")
	open_out.write(taxa_title + "\tsource\tmsp_name\tmsp_taxa_name\tmsp_taxa_id\n")
	for myid in sorted(genes.keys()):
		info = genes[myid].split("\t")
		mystr = genes[myid]
		map_type = info[titles["map_type"]]
		if map_type == "UniRef90_uncharacterized" or map_type == "UniRef90_characterized": 	# UniRef90_strong_homology
			# save the taxonomy of the best UniRef90 hit for each protein as a putative taxon only. These can be used for defining the LCA taxon for an MSP, but we don't trust them much on a per-protein basis UNLESS the hit was a strong hit.
			if myid in msp: # grouped into MSP
				mymsp, myinfo = msp[myid].split("\n")
				msp_taxa = "NA"
				msp_taxa_id = "NA"
				if mymsp in anns:
					tmp1 = anns[mymsp].split("\t")
					msp_taxa = tmp1[2]
					msp_taxa_id = tmp1[1]
				mystr = genes[myid] + "\tUniRef90\t" + mymsp + "\t" + msp_taxa + "\t" + msp_taxa_id
				if info[titles["taxa_rank"]] == "Unclassified" or re.search("Kingdom", info[titles["taxa_rank"]]):
					if mymsp in anns:
						tmp1 = anns[mymsp].split("\t")
						info[titles["taxa_id"]] = tmp1[1]
						info[titles["taxa_name"]] = tmp1[2]
						info[titles["taxa_rank"]] = tmp1[3]
						info[titles["taxa_lineage"]] = tmp1[4]
						tmp_line = "\t".join(info)
						mystr = tmp_line + "\tUniRef90_unclassified_MSP\t" + mymsp + "\t" + tmp1[2] + "\t" + tmp1[1]
					else:
						# debug
						print("No MSP annotation info!\t" + mymsp)
						if info[titles["taxa_rank"]] == "Unclassified" or re.search("Kingdom", info[titles["taxa_rank"]]):	
							info[titles["taxa_id"]] = "NA"
							info[titles["taxa_name"]] = "NA"
							info[titles["taxa_rank"]] = "Unclassified"
							info[titles["taxa_lineage"]] = "NA"
							tmp_line = "\t".join(info)
							mystr = tmp_line + "\tUniRef90_unclassified_non-MSP\t" + "msp_unknown\tNA\tNA"
			else:	# not grouped into MSP
				mystr = mystr + "\tUniRef90\t" + "msp_unknown\tNA\tNA"
				if info[titles["taxa_rank"]] == "Unclassified" or re.search("Kingdom", info[titles["taxa_rank"]]):
					info[titles["taxa_id"]] = "NA"
					info[titles["taxa_name"]] = "NA"
					info[titles["taxa_rank"]] = "Unclassified"
					info[titles["taxa_lineage"]] = "NA"
					tmp_line = "\t".join(info)
					mystr = tmp_line + "\tUniRef90_unclassified_non-MSP\t" + "msp_unknown\tNA\tNA"
		else: # UniRef90_unknown
			# when multiple weak votes align inside of an MSP, we can trust them more. We should not trust the taxonomic assignment to an isolate weak homolog very much
			if myid in msp: # grouped into MSP
				mymsp, myinfo = msp[myid].split("\n")
				if mymsp in anns:
					tmp1 = anns[mymsp].split("\t")
					info[titles["taxa_id"]] = tmp1[1]
					info[titles["taxa_name"]] = tmp1[2]
					info[titles["taxa_rank"]] = tmp1[3]
					info[titles["taxa_lineage"]] = tmp1[4]
					tmp_line = "\t".join(info)
					mystr = tmp_line + "\tUniRef90_unclassified_MSP\t" + mymsp + "\t" + tmp1[2] + "\t" + tmp1[1]
				else:
					# debug
					print("No MSP annotation info!\t" + mymsp)
					info[titles["taxa_id"]] = "NA"
					info[titles["taxa_name"]] = "NA"
					info[titles["taxa_rank"]] = "Unclassified"
					info[titles["taxa_lineage"]] = "NA"
					tmp_line = "\t".join(info)
					mystr = tmp_line + "\tUniRef90_unclassified_non-MSP\t" + "msp_unknown\tNA\tNA"
			else:	# not grouped into MSP
				info[titles["taxa_id"]] = "NA"
				info[titles["taxa_name"]] = "NA"
				info[titles["taxa_rank"]] = "Unclassified"
				info[titles["taxa_lineage"]] = "NA"
				tmp_line = "\t".join(info)
				mystr = tmp_line + "\tUniRef90_unclassified_non-MSP\t" + "msp_unknown\tNA\tNA"
		open_out.write(mystr + "\n")
	# foreach gene
	open_out.close()
# func: assign_msp


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start extract_MSPminer_annotation.py -a " + values.annotation + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	anns = collect_taxonomy_info (values.annotation, values.type)
	msp = collect_msp_info (values.msp)
	assign_msp (values.gene, msp, anns, values.output)
	sys.stderr.write("Get info ......done\n")


	sys.stderr.write("### Finish extract_MSPminer_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

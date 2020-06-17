#!/usr/bin/env python

"""
MetaWIBELE: ddi_DOMINE_SIFTS module
Based on SIFTS DBs, divide DDIs into known human+bacterial DDIs vs. interologs 

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
Based on SIFTS DBs, divide DDIs into known human+bacterial DDIs vs. interologs
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--ddi-annotation",
	                    help='input the annotated DDI file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output summary detail file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect info
#==============================================================
def collect_taxa_list (taxa_file): # uniprot_taxaID_bac-arc-vir.tsv
	taxa = {}
	for line in utils.gzip_bzip2_biom_open_readlines (taxa_file):
		line = line.strip()
		if not len(line):
			continue
		taxa[line] = ""
	# foreach line
	
	return taxa
# collect_taxa_list


def collect_pfam_info (pfam_file):	# pdb_chain_pfam.tsv 
	pfams = {}
	for line in utils.gzip_bzip2_biom_open_readlines (pfam_file):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line) or re.search("^PDB", line):
			continue
		info = line.split("\t")
		mypdb = info[0]
		mychain = info[1]
		mypfam = info[-2]
		if not mypdb in pfams:
			pfams[mypdb] = {}
		if not mychain in pfams[mypdb]:
			pfams[mypdb][mychain] = {}
		pfams[mypdb][mychain][mypfam] = ""
	# foreach line
	
	return pfams
# function collect_pfam_info


def collect_taxanomy_info (map_file):	# pdb_chain_taxonomy.tsv 
	taxa = {}
	taxa_hit = {}
	for line in utils.gzip_bzip2_biom_open_readlines (map_file): 
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line) or re.search("^PDB", line):
			continue
		info = line.split("\t")
		mypdb = info[0]
		mychain = info[1]
		tax = info[2]
		if not mypdb in taxa:
			taxa[mypdb] = {}
		if not mypdb in taxa_hit:
			taxa_hit[mypdb] = {}
		taxa_hit[mypdb][tax] = ""
		if not mychain in taxa[mypdb]:
			taxa[mypdb][mychain] = {}
		taxa[mypdb][mychain][tax] = ""
	# foreach line

	#taxa_flt = {}
	#for mypdb in taxa:
	#	if mypdb in taxa_hit:
	#		if filter_flag == "yes":
	#			if not "9606" in taxa_hit[mypdb]:	# human PDB
	#				continue
	#		# filtering human PDB
	#		if not mypdb in taxa_flt:
	#			taxa_flt[mypdb] = {}
	#		for mychain in taxa[mypdb]:
	#			if not mychain in taxa_flt[mypdb]:
	#				taxa_flt[mypdb][mychain] = {}
	#			for mytax in taxa[mypdb][mychain]:
	#				taxa_flt[mypdb][mychain][mytax] = ""
	#		# foreach chain
	#	# if hit 
	# foreach PDB
	return taxa
# function collect_taxonomy_info


def collect_interaction_info (ann_file):	# summary_DOMINE_peptide.denovo.ann.tsv
	interact = {}
	titles = {}
	title = ""
	open_file = open(ann_file, "r")
	line = open_file.readline()
	line = line.strip()
	title = line
	info = line.split("\t")
	myindex = 0
	myflag = info[0]
	while myindex < len(info):
		titles[info[myindex]] = myindex
		myindex = myindex + 1
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid1 = info[0]
		pfam1 = info[titles["Pfam1_ID"]]
		pfam2 = info[titles["Pfam2_ID"]]
		myid2 = pfam1 + "\t" + pfam2
		if not myid1 in interact:
			interact[myid1] = {}
		interact[myid1][myid2] = line
	# foreach line
	open_file.close()
	return interact, title, myflag
# collect_interaction_info


#==============================================================
# assign structure-supported interactions to DDIs 
#==============================================================
def assign_interaction (filter_flag, taxa1, taxa2, cutoff, pfams, pdb_taxa, interact, title, id_flag, outfile):
	# SIFTS structure supporting
	PDB = {}
	PDB_all = {}
	for mypdb in pdb_taxa:
		chain_num = len(pdb_taxa[mypdb].keys())
		if chain_num < 2:	# only one chain
			continue
		for mychain1 in sorted(pdb_taxa[mypdb].keys()):
			mytaxa1 = ";".join(sorted(pdb_taxa[mypdb][mychain1].keys()))
			for mychain2 in sorted(pdb_taxa[mypdb].keys()):
				if mychain1 == mychain2:
					continue
				mytaxa2 = ";".join(sorted(pdb_taxa[mypdb][mychain2].keys()))
				taxa_pair = mytaxa1 + "\t" + mytaxa2
				chain_pair = mychain1 + ":" + mychain2
				# check specified chain structure
				support = 0
				if filter_flag == "yes":
					flag1 = 0
					flag2 = 0
					for taxa_1 in pdb_taxa[mypdb][mychain1].keys():
						if taxa_1 in taxa1:
							flag1 = 1
							break
					for taxa_2 in pdb_taxa[mypdb][mychain2].keys():
						if taxa_2 in taxa2:
							flag2 = 1
							break
					if flag1 * flag2 == 0:
						flag1 = 0
						flag2 = 0
						for taxa_1 in pdb_taxa[mypdb][mychain2].keys():
							if taxa_1 in taxa1:
								flag1 = 1
								break
						for taxa_2 in pdb_taxa[mypdb][mychain1].keys():
							if taxa_2 in taxa2:
								flag2 = 1
								break
					if flag1 * flag2 == 1:	# have specified PDB structure
						support = 1
				else:
					if mytaxa1 == mytaxa2 and not re.search(";", mytaxa1):	# at leat two species intereaction
						support = 0
					else:
						support = 1
				# check pfam
				if mypdb in pfams:
					if mychain1 in pfams[mypdb] and mychain2 in pfams[mypdb]:
						for mypfam1 in sorted(pfams[mypdb][mychain1].keys()):
							for mypfam2 in sorted(pfams[mypdb][mychain2].keys()):
								pfam_pair = mypfam1 + "\t" + mypfam2
								if not pfam_pair in PDB_all:
									PDB_all[pfam_pair] = {}
								PDB_all[pfam_pair][taxa_pair] = mypdb + "\t" + chain_pair
								if support == 1: # supporting by SIFTS
									if mytaxa1 == mytaxa2 and not re.search(";", mytaxa1):  # at leat two species intereaction								
										continue
									if not pfam_pair in PDB:
										PDB[pfam_pair] = {}
									PDB[pfam_pair][taxa_pair] = mypdb + "\t" + chain_pair
							# foreach pfam2
						# foreach pfma1
					# if chain1 and chain2
				# if pfam
			# foreach chain2
		# foreach chain1
	# foreach PDB

	# check interaction
	number = {}
	anns = {}
	for pepid in sorted(interact.keys()):
		flag = 0
		for pfam_pair in sorted(interact[pepid].keys()):
			if pfam_pair in PDB_all:
				flag = 1
				if not pfam_pair in PDB:
					myinter = "SIFTS_not_supporting"
					if cutoff == "loose":
						for mytaxa in sorted(PDB_all[pfam_pair].keys()):
							if not pepid in anns:
								anns[pepid] = {}
							if not pfam_pair in anns[pepid]:
								anns[pepid][pfam_pair] = []
							anns[pepid][pfam_pair].append(interact[pepid][pfam_pair] + "\t" + PDB_all[pfam_pair][mytaxa] + "\t" + mytaxa + "\t" + myinter)
				else:
					myinter = "SIFTS_supporting"
					for mytaxa in sorted(PDB[pfam_pair].keys()):
						if not pepid in anns:
							anns[pepid] = {}
						if not pfam_pair in anns[pepid]:
							anns[pepid][pfam_pair] = []
						anns[pepid][pfam_pair].append(interact[pepid][pfam_pair] + "\t" + PDB[pfam_pair][mytaxa] + "\t" + mytaxa + "\t" + myinter)
				if not myinter in number:
					number[myinter] = {}
				number[myinter][pepid] = pfam_pair
				# foreach taxa pair
			# if pfam_pair in SIFTS
		# foreach pair
		if flag == 0: # no SIFTS hits
			myinter = "SIFTS_unknown"
			if not myinter in number:
				number[myinter] = {}
			number[myinter][pepid] = pfam_pair
	# foreach cluster

	outs = {}
	details = {}
	outfile1 = re.sub(".detail.tsv", ".tsv", outfile)
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile, "w")
	open_out1.write(title + "\t" + "PDB\t" + "Chain_pair\t" + "Tax1_ID\tTax2_ID\tSIFTS_type" + "\n")
	open_out2.write(id_flag + "\ttype\tdetail\tPDB\n")
	for mypep in sorted(anns.keys()):
		for mypair in sorted(anns[mypep].keys()):
			tmp = re.sub("\t", ":", mypair)
			for item in anns[mypep][mypair]:
				info = item.split("\t")
				mytype = info[-1]
				myann = info[1] + "\t" + info[2] + "\t" + info[3] + "\t" + info[4]
				if not myann in outs:
					outs[myann] = {}
				outs[myann][mypep] = ""
				if not mypep in details:
					details[mypep] = {}
				details[mypep][tmp] = info[-5] + "," + info[-4]
				open_out1.write(item + "\n")
			# foreach taxa_pair
		# foreach pfam_pair
	# foreah pepid
	open_out1.close()
	for mypep in sorted(details.keys()):
		ddi_str = ""
		pdb_str = ""
		for myddi in sorted(details[mypep].keys()):
			ddi_str = ddi_str + myddi + ";"
			pdb_str = pdb_str + details[mypep][myddi] + ";"
		# foreach DDI
		ddi_str = re.sub(";$", "", ddi_str)
		pdb_str = re.sub(";$", "", pdb_str)
		open_out2.write(mypep + "\tSIFTS_interaction\t" + ddi_str + "\t" + pdb_str + "\n")
	# foreach cluster
	open_out2.close()

	outfile3 = re.sub(".tsv", ".info.tsv", outfile1)
	open_out3 = open(outfile3, "w")
	open_out3.write("Pfam1_ID\tPfam1_ann\tPfam2_ID\tPfam2_ann\tnumber\n")
	orders = {}
	for myid in outs.keys():
		mynum = len(outs[myid].keys())
		if not mynum in orders:
			orders[mynum] = {}
		orders[mynum][myid] = ""
	# foreach pair
	for mynum in sorted(orders.keys(), key=int, reverse=True):
		for myid in sorted(orders[mynum].keys()):
			open_out3.write(myid + "\t" + str(mynum) + "\n")
		# foreach id
	# foreach number
	open_out3.close()

	"""
	outfile2 = re.sub(".tsv", ".plot.tsv", outfile1)
	open_out = open(outfile2, "w")
	open_out.write("type\t" + id_flag + "\tPfam1_ID\tPfam2_ID\n")
	for mytype in sorted(number.keys()):
		for mypep in sorted(number[mytype].keys()):
			open_out.write(mytype + "\t" + mypep + "\t" + number[mytype][mypep] + "\n")
	open_out.close()
	"""

# assign_annotation


#==============================================================
###########  Main processing ############
#=============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start ddi_DOMINE_SIFTS.py -i " + values.ddi_annotation + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	taxa1 = collect_taxa_list (config.microbiome_taxa)
	taxa2 = collect_taxa_list (config.mammalia_taxa)
	pfams = collect_pfam_info (config.pdb_pfam)
	pdb_taxa = collect_taxanomy_info (config.pdb_taxonomy)
	interact, title, id_flag = collect_interaction_info (values.ddi_annotation)
	sys.stderr.write("Get info ......done\n")
	
	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign SIFTS structure to peptide families ......starting\n")
	cutoff = "stringent"
	filter = "yes"
	assign_interaction ("yes", taxa1, taxa2, "stringent", pfams, pdb_taxa, interact, title, id_flag, values.output)
	sys.stderr.write("\nAssign SIFTS structure to peptide families ......done\n")

	sys.stderr.write("### Finish ddi_DOMINE_SIFTS.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

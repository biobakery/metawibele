#!/usr/bin/env python

"""
MetaWIBELE: mspminer_msp_taxonomy_annotation module
Taxonomy annotation for MSPs 

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
from collections import Counter

try:
	from metawibele import config
	from metawibele.common import utils
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Taxonomy annotation for MSPs 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--msp-annotation",
	                    help='input the MSPs annotation file',
	                    required=True)
	parser.add_argument('-t', "--homology-type",
	                    help='specify the type of uniref90 homology',
	                    choices=["no", "UniRef90_strong_homology", "UniRef90_homology"],
	                    required=True,
	                    default="UniRef90_homology")
	parser.add_argument('-o', "--output",
	                    help='output taxonomy annotation of MSP file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# extract taxon for MSPs
#==============================================================
def collect_taxonomy_info (map_file):  
	taxa_map = {}
	titles = {}
	for line in utils.gzip_bzip2_biom_open_readlines (map_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^Taxon", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		# title line
		taxa_id = info[titles["Taxon"]]
		taxa_name = info[titles["Scientific_name"]]
		taxa_rank = info[titles["Rank"]]
		taxa_lineage = info[titles["Lineage"]]
		myid = taxa_name
		tmp = taxa_lineage.split("|")
		if re.search("__", tmp[-1]):
			mym = re.search("([^\_]+)__([\S]+)", tmp[-1])
			myrank = mym.group(1)
			if myrank == "k":
				myrank = "Kingdom"
			if myrank == "p":
				myrank = "Phylum"
			if myrank == "c":
				myrank = "Class"
			if myrank == "o":
				myrank = "Order"
			if myrank == "f":
				myrank = "Family"
			if myrank == "g":
				myrank = "Genus"
			if myrank == "s":
				myrank = "Species"
			if myrank == "t":
				myrank = "Terminal"
			taxa_rank = myrank
			myid = mym.group(2)
		taxa_map[myid] = taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage + "\t" + taxa_rank + "__" + taxa_name
	# foreach line
	
	return taxa_map
# collect_taxonomy_info



#==============================================================
# extract taxon for MSPs
#==============================================================
def extract_taxon_info (msp_file, homology):  
	# collect info
	msp = {}
	taxa = {}
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
		msp_name = info[titles["msp_name"]]
		myclass = info[titles["class"]]
		module_name = info[titles["module_name"]]
		gene_name = info[titles["gene_name"]]
		taxID = info[titles["taxa_id"]]
		taxa_name = info[titles["taxa_name"]]
		taxa_rank = info[titles["taxa_rank"]]
		taxa_lineage = info[titles["taxa_lineage"]]
		#organism = info[titles["organism"]]
		map_type = info[titles["map_type"]]
		identity = info[titles["identity"]]
		coverage = info[titles["mutual_coverage"]]
		if homology != "no":
			if homology == "UniRef90_strong_homology":
				if map_type != "UniRef90_characterized" and map_type != "UniRef90_uncharacterized":
					taxID = "NA"
					taxa_name = "NA"
					taxa_rank = "Unclassified"
					taxa_lineage = "NA"
					#organism = "NA"
			if homology == "UniRef90_homology":
				if map_type == "UniRef90_worse_homology":
					taxID = "NA"
					taxa_name = "NA"
					taxa_rank = "Unclassified"
					taxa_lineage = "NA"
					#organism = "NA"
				else:
					if not(float(identity) >= 50 and float(coverage) >= 0.8):	# UniRef50-like protein
						taxID = "NA"
						taxa_name = "NA"
						taxa_rank = "Unclassified"
						taxa_lineage = "NA"
						#organism = "NA"
		if not msp_name in msp:
			msp[msp_name] = {}
		if not myclass in msp[msp_name]:
			msp[msp_name][myclass] = {}
		msp[msp_name][myclass][gene_name] = taxa_lineage
		if not re.search("shared_accessory", myclass):
			if not "all" in msp[msp_name]:
				msp[msp_name]["all"] = {}
			msp[msp_name]["all"][gene_name] = taxa_lineage
		if taxa_lineage != "NA":
			tmp = taxa_lineage.split("|")
			for item in tmp:
				if re.search("__", item):
					mym = re.search("^([^\_]+)__([\S]+)", item)
					myrank = mym.group(1)
					myname = mym.group(2)
					if myname == "NA":
						continue
					if myrank == "k":
						myrank = "Kingdom"
					if myrank == "p":
						myrank = "Phylum"
					if myrank == "c":
						myrank = "Class"
					if myrank == "o":
						myrank = "Order"
					if myrank == "f":
						myrank = "Family"
					if myrank == "g":
						myrank = "Genus"
					if myrank == "s":
						myrank = "Species"
					if myrank == "t":
						myrank = "Terminal"
					#myname = re.sub("_", " ", myname)
					if not myrank in taxa:
						taxa[myrank] = {}
					#if not gene_name in taxa[myrank]:
					taxa[myrank][gene_name] = myname
				# if taxon
			# foreach taxon rank
		# if taxon exists
	# foreach gene
	open_file.close()
	return msp, taxa
# extract_taxon_info


#==============================================================
# taxonomy annotation for MSPs
#==============================================================
def taxonomy_annotation (msp, taxa, normalized_known, min_known, min_cutoff, cutoff_type, taxa_map, outfile):
	open_out = open(outfile, "w")
	open_out.write("msp_name\tclass\tnum_genes\tnum_unclassified_genes\tfraction_unclassified_genes\tnum_hit_genes\tfraction_hit_genes\ttaxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\ttaxa_rank_name\n")
	#taxa_level = ["Terminal", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]
	taxa_level = ["Terminal", "Species", "Genus", "Family", "Order", "Class", "Phylum"]
	for msp_name in sorted(msp.keys()):
		for myclass in sorted(msp[msp_name].keys()):
			taxa_flag = 0
			taxa_name = "NA"
			hit_num = "NA"
			hit_per = "NA"
			total = 0
			unknown = 0
			known = 0
			taxa_num = {}
			for mygene in msp[msp_name][myclass].keys():
				myline = msp[msp_name][myclass][mygene]
				total = total + 1
				if myline == "NA":
					unknown = unknown + 1
				for mytaxa in taxa.keys():
					if mygene in taxa[mytaxa]:
						myname = taxa[mytaxa][mygene]
						if not mytaxa in taxa_num:
							taxa_num[mytaxa] = []
						taxa_num[mytaxa].append(myname)
				# foreach taxon rank
			# foreach gene
			known = total - unknown
			if float(known)/float(total) < float(min_known):
				# very few classified genes in the MSP and fail in assigning taxon to whole group
				taxa_flag = 1
				taxa_name = "NA"
				hit_num = 0
				hit_per = 0
			# unclassified
			for mylevel in taxa_level:
				if taxa_flag == 1:
					break
				if mylevel in taxa_num:
					tmp = Counter(taxa_num[mylevel])
					first_tax = "NA"
					first_num = 0
					second_tax = "NA"
					second_num = 0
					try:
						for (mytax, tax_num) in sorted(tmp.iteritems(), key = lambda d:d[1], reverse = True): # python2
							if first_tax == "NA": # most dominant
								first_tax = mytax
								first_num = int(tax_num)
								continue
							if first_tax != "NA" and second_tax == "NA": # second dominant
								second_tax = mytax
								second_num = int(tax_num)
								break
					except:
						for (mytax, tax_num) in sorted(tmp.items(), key = lambda d:d[1], reverse = True): # python3
							if first_tax == "NA": # most dominant
								first_tax = mytax
								first_num = int(tax_num)
								continue
							if first_tax != "NA" and second_tax == "NA": # second dominant
								second_tax = mytax
								second_num = int(tax_num)
								break
					# foreach type of taxon
					if cutoff_type == "most":
						mydiff = first_num
					if cutoff_type == "diff":
						mydiff = first_num - second_num
					if normalized_known == "yes":	# normalized to known taxon level
						myper = float(mydiff) / float(known)
						myper1 = float(first_num) / float(known)
					else:
						myper = float(mydiff) / float(total)
						myper1 = float(first_num) / float(total)
					if myper >= float(min_cutoff):	# above the consistency cutoff
						taxa_flag = 1
						taxa_name = first_tax
						hit_num = first_num
						hit_per = myper1
				# if taxa_num
			# foreach taxon level
			unknown_per = float(unknown) / float(total)
			if taxa_flag == 0:	# still unclassified
				taxa_name = "NA"
				hit_num = first_num
				hit_per = myper1
				# debug
				print("No taxon that is sig. more dominant than the second one!\t" + msp_name + "\t" + myclass + "\t" + first_tax + "\t" + str(first_num) + "\t" + second_tax + "\t" + str(second_num))
				if float(myper1) > 0.50:
					print("Assign the MSP using the most dominant taxon given most members are classified\t" + msp_name + "\t" + myclass + "\t" + first_tax + "\t" + str(myper1))
					taxa_name = first_tax
			taxa_info = "NA\tNA\tUnclassified\tNA\tNA"
			if taxa_name in taxa_map:
				taxa_info = taxa_map[taxa_name]
				tmp = taxa_info.split("\t")
				if re.search("Superkingdom", tmp[2]) or re.search("Kingdom", tmp[2]): # kindom as unclassified
					taxa_info = "NA\tNA\tUnclassified\tNA\tNA"
			mystr = msp_name + "\t" + myclass + "\t" + str(total) + "\t" + str(unknown) + "\t" + str(unknown_per) + "\t" + str(hit_num) + "\t" + str(hit_per) + "\t" + taxa_info
			open_out.write(mystr + "\n")
		# foreach class of gene
	# foreach msp
	open_out.close()
# taxonomy_annotation 


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start mspminer_msp_taxonomy_annotation.py -a " + values.msp_annotation + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	taxa_map = collect_taxonomy_info (config.taxonomy_database)
	msp, taxa = extract_taxon_info (values.msp_annotation, values.homology_type)
	normalized_known = "yes"
	cutoff_type = "diff"
	taxonomy_annotation (msp, taxa, normalized_known, config.tshld_unclassified, config.tshld_diff, cutoff_type, taxa_map, values.output)
	sys.stderr.write("Get info ......done\n")

	sys.stderr.write("### Finish mspminer_msp_taxonomy_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

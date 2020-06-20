#!/usr/bin/env python
##########################################################################
# Function: Format the taxonomy lineage information from UniProt taxonomy file 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 10/01/2018
##########################################################################
import sys
import os
import re
import argparse


try:
	from metawibele.common import utils
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
            " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Format the taxonomy lineage information from UniProt taxonomy file 
"""

def get_args ():
	parser=argparse.ArgumentParser()
	parser.add_argument('-i', "--input",
						help='the input UniProt taxonomy file', 
						required=True)
	parser.add_argument('-o', "--output",
						help='the path of output folder', 
						required=True)
	values=parser.parse_args()
	return values
# get_args


#==============================================================
# Format taxonomy info
#==============================================================
def format_taxonomy_info (taxafile, output_path):	# uniprot_taxonomy.tsv
	titles = {}
	taxa_info = {}
	kingdoms = {}
	phylums = {}
	classes = {}
	orders = {}
	families = {}
	gena = {}
	species = {}
	kingdom = ["Superkingdom", "Kingdom", "Subkingdom"]
	phylum = ["Superphylum", "Phylum", "Subphylum"]
	clas = ["Superclass", "Class", "Subclass", "Infraclass", "Cohort"]
	order = ["Superorder", "Order", "Suborder", "Infraorder", "Parvorder"]
	family = ["Superfamily", "Family", "Subfamily", "Tribe", "Subtribe"]
	genus = ["Genus", "Subgenus"]
	specie = ["Species", "Species group", "Species subgroup", "Species subsgroup", "Subspecies", "Forma", "Varietas"]
		
	# collect info
	title_num = 0
	for line in utils.gzip_bzip2_biom_open_readlines (taxafile): 
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		info = line.split("\t")
		if re.search("^Taxon", line):
			myindex = 0
			while myindex < len(info):
				item = info[myindex]
				titles[item] = myindex
				myindex = myindex + 1
			# foreach item
			title_num = len(info)
			continue	
		
		# debug
		if len(info) < title_num:
			mynum_tmp = len(info)
			while mynum_tmp < title_num:
				info.append("")
				mynum_tmp = mynum_tmp + 1		
		#print(str(len(info)) + "\t" + line)

		mytaxa = info[titles["Taxon"]]
		myname = info[titles["Scientific name"]]
		if myname == "":
			continue
		myrank = info[titles["Rank"]]
		myline = info[titles["Lineage"]]
		mypar = info[titles["Parent"]]
		mylevel = re.sub("\s+", "_", myname)
		if myrank in kingdom:
			mylevel = "k__" + re.sub("\s+", "_", myname)
			kingdoms[myname] = mytaxa
		if myrank in phylum:
			mylevel = "p__" + re.sub("\s+", "_", myname)
			phylums[myname] = mytaxa
		if myrank in clas:
			mylevel = "c__" + re.sub("\s+", "_", myname)
			classes[myname] = mytaxa
		if myrank in order:
			mylevel = "o__" + re.sub("\s+", "_", myname)
			orders[myname] = mytaxa
		if myrank in family:
			mylevel = "f__" + re.sub("\s+", "_", myname)
			families[myname] = mytaxa
		if myrank in genus:
			mylevel = "g__" + re.sub("\s+", "_", myname)
			gena[myname] = mytaxa
		if myrank in specie:
			mylevel = "s__" + re.sub("\s+", "_", myname)
			species[myname] = mytaxa
		taxa_info[mytaxa] = mytaxa + "\t" + myname + "\t" + myrank + "\t" + myline + "\t" + mylevel + "\t" + mypar	
	# foreach line
	
	## format taxonomy info
	mic_types = ["k__Viruses", "k__Bacteria", "k__Archaea", "k__Fungi"]
	mammalia_types = ["c__Mammalia"]
	mic_ids = {}
	mammalia_ids = {}
	outfile = os.path.join(output_path, "uniprot_taxonomy.tsv")
	open_out = open(outfile, "w")
	open_out.write("Taxon\tScientific_name\tRank\tLineage\tParent\n")
	for mytaxa in sorted(taxa_info.keys()):
		mytaxa, myname, myrank, myline, mylevel, mypar = taxa_info[mytaxa].split("\t")
		if myrank == "":
			mylevel = "t__" + mylevel
		info = myline.split("; ")
		mystr = ""
		if len(info) < 1:
			mystr = mylevel
		else:
			mystr = ""
			for item in info:
				if item in kingdoms:
					item = "k__" + item
				if item in phylums:
					item = "p__" + item
				if item in classes:
					item = "c__" + item
				if item in orders:
					item = "o__" + item
				if item in families:
					item = "f__" + item
				if item in gena:
					item = "g__" + item
				if item in species:
					item = "s__" + item
				if item != "":
					item = re.sub("\s+", "_", item)
					mystr = mystr + item + "|"
			# foreach item
			# check parent
			tmp_info = mystr.split("|")
			if mypar in taxa_info:
				mypar_level = taxa_info[mypar].split("\t")[-2]
				if not mypar_level in tmp_info:	# skipped the nearest ancester
					mystr = mystr + mypar_level + "|"
			mystr = mystr + mylevel
		# else
		# rename empoty myrank
		if myrank == "":
			if re.search("t__", mylevel):
				myrank = "Terminal"
			if re.search("s__", mylevel):
				myrank = "Species"
			if re.search("g__", mylevel):
				myrank = "Genus"
			if re.search("f__", mylevel):
				myrank = "Family"
			if re.search("o__", mylevel):
				myrank = "Order"
			if re.search("c__", mylevel):
				myrank = "Class"
			if re.search("p__", mylevel):
				myrank = "Phylum"
			if re.search("k__", mylevel):
				myrank = "Kingdom"
		open_out.write(mytaxa + "\t" + myname + "\t" + myrank + "\t" + mystr + "\t" + mypar + "\n")

		for i in mic_types:
			if re.search(i, mystr):
				mic_ids[mytaxa] = ""
		for i in mammalia_types:
			if re.search(i, mystr):
				mammalia_ids[mytaxa] = ""
	# foreach taxon
	open_out.close()
	os.system("gzip " + outfile)

	# output microbiome and mammalia taxa ids
	outfile = os.path.join(output_path, "uniprot_taxaID_microbiome.txt")
	open_out = open(outfile, "w")
	for myid in sorted(mic_ids.keys()):
		open_out.write(myid + "\n")
	open_out.close()
	os.system("gzip " + outfile)
	outfile = os.path.join(output_path, "uniprot_taxaID_mammalia.txt")
	open_out = open(outfile, "w")
	for myid in sorted(mammalia_ids.keys()):
		open_out.write(myid + "\n")
	open_out.close()
	os.system("gzip " + outfile)

# format_taxanomy_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():

	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start prepare_uniprot_taxonomy.py -i " + values.input + " ####\n")
	
	### Extraction ###
	sys.stderr.write("Extract taxonomic info......starting\n")
	format_taxonomy_info (values.input, values.output)
	sys.stderr.write("Extract taxonomic info......done\n")
	
	sys.stderr.write("### Finish prepare_uniprot_taxonomy.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

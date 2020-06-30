#!/usr/bin/env python

"""
MetaWIBELE: ddi_DOMINE_protein module
Summary the domain-domian interaction based on DOMINE database 

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
Summary the domain-domian interaction based on DOMINE database
"""

def get_args ():
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
	                    default="interpro.PfamDomain.tsv")
	parser.add_argument('-p', "--path",
	                    help='input the path of annotation file',
	                    required=True)
	parser.add_argument('-f', "--filter",
	                    help='whether filter for human-microbiome DDIs',
	                    choices=["yes", "no"],
	                    required=True,
	                    default="no")
	parser.add_argument('-s', "--suffix",
	                    help='specify the name suffix of DDI output file',
	                    required=True,
	                    default="DDI.tsv")

	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect info
#==============================================================
def collect_pfam_ann (pfamfile):   # Pfam_ann.tsv 
	pfams = {}
	for line in utils.gzip_bzip2_biom_open_readlines (pfamfile):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		info = line.split("\t")
		pfam = info[0]
		ann = info[1]
		if not pfam in pfams:
			pfams[pfam] = ann
    # foreach line
	
	return pfams
#collect_pfam_ann


def collect_pfam_info (pfam_file):	# split1.interpro.PfamDomain.tsv
	peptide = {}
	titles = {}
	open_file = open(pfam_file, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_ID , line):
			for item in info:
				titles[item] = info.index(item)
			continue
		info = line.split("\t")
		myid = info[titles[utilities.PROTEIN_ID]]
		pfam = info[titles["Pfam"]]
		ann = info[titles["Description"]]
		if pfam == "NA":
			continue
		if not myid in peptide:
			peptide[myid] = {}
		peptide[myid][pfam] = ann
	# foreach line
	open_file.close()
	
	return peptide
# function collect_pfam_info


def collect_mapping_info (map_file):	# uniprot_human_pfam.tsv
	maps = {}
	titles = {}
	for line in utils.gzip_bzip2_biom_open_readlines (map_file):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		if re.search("^Pfam\t", line):
			info = line.split("\t")
			myindex = 0
			while myindex < len(info):
				titles[info[myindex]] = myindex
				myindex = myindex + 1
			continue
		info = line.split("\t")
		pfam = info[titles["Pfam"]]
		organism = info[titles["Organism"]]
		taxaID = info[titles["NCBI_TaxID"]]
		if not pfam in maps:
			maps[pfam] = {}
		maps[pfam][taxaID] = ""
	# foreach line

	# human pfams
	human_pfam = {}
	for mypfam in maps:
		if "9606" in maps[mypfam]:
			human_pfam[mypfam] = ";".join(sorted(maps[mypfam].keys()))
	return  human_pfam
# function collect_mapping_info

def collect_interaction_info (int_file, filter_flag, human_pfam):	# INTERACTION.txt
	interact = {}
	for line in utils.gzip_bzip2_biom_open_readlines (int_file): 
		line = line.strip()
		if not len(line):
			continue
		info = line.split("|")
		id1 = info[0]
		id2 = info[1]
		level = info[-2]
		#if level == "NA":
		#	continue
		if id1 == id2:	# self-interaction
			continue
		if filter_flag == "yes":
			if not id1 in human_pfam and not id2 in human_pfam: # no interaction with human pfam
				continue
		if not id1 in interact:
			interact[id1] = {}
		interact[id1][id2] = level
		if not id2 in interact:
			interact[id2] = {}
		interact[id2][id1] = level
	# foreach line
	
	return interact
# collect_interaction_info


#==============================================================
# assign domain interaction 
#==============================================================
def assign_interaction (filter_flag, spe_level, pfams, human_pfam, interact, peptide, outfile, outfile_detail):
	outs = {}
	outs_info = {}
	anns = {}
	for pepid in sorted(peptide.keys()):
		for pfam1 in sorted(peptide[pepid].keys()):
			if not pfam1 in interact:
				continue
			for pfam2 in sorted(interact[pfam1].keys()):
				mylevel = interact[pfam1][pfam2]
				if filter_flag == "yes":
					if not pfam2 in human_pfam:	# interacted domain not in human
						continue
				if not pepid in outs:
					outs[pepid] = {}
				mytype = "DOMINE_interaction"
				myid_tmp = mytype + "\t" + mylevel
				if not myid_tmp in outs[pepid]:
					outs[pepid][myid_tmp] = pfam1 + ":" + pfam2
				else:
					outs[pepid][myid_tmp] = outs[pepid][myid_tmp] + ";" + pfam1 + ":" + pfam2
				item = mytype + "\t" + mylevel + "\t" + pfam1 + "\t" + pfam2
				item_v = mytype + "\t" + mylevel + "\t" + pfam2 + "\t" + pfam1
				if not pepid in outs_info:
					outs_info[pepid] = {}
				if not item_v in outs_info[pepid]:
					outs_info[pepid][item] = ""
				if not pepid in anns:
					anns[pepid] = {}
				anns[pepid][mylevel] = ""
			# foreach pfam2
		# foreach pfam1
	# foreach cluster

	open_out = open(outfile, "w")
	open_out.write(utilities.PROTEIN_ID + "\tType\tInteraction\tPfam1_ID\tPfam2_ID\tPfam1_ann\tPfam2_ann\n")
	for mypep in sorted(outs_info.keys()):
		for myinfo in sorted(outs_info[mypep].keys()):
			tmp = myinfo.split("\t")
			myann1 = "NA"
			myann2 = "NA"
			if tmp[-2] in pfams:
				myann1 = pfams[tmp[-2]]
			if tmp[-1] in pfams:
				myann2 = pfams[tmp[-1]]
			open_out.write(mypep + "\t" + myinfo + "\t" + myann1 + "\t" + myann2 + "\n")
		# foreach item
	# foreah pepid
	open_out.close()

	open_out = open(outfile_detail, "w")
	open_out.write(utilities.PROTEIN_ID + "\ttype\tdetail\tannotation\tinteraction\n")
	for mypep in sorted(outs.keys()):
		pfam_info = ""
		ann_info = ""
		level_info = ""
		for myid in sorted(outs[mypep].keys()):
			mytype, mylevel  = myid.split("\t")
			if spe_level != "no":
				if mylevel != "NA":
					if mylevel != spe_level:	# not specified interaction level
						continue
			myann = "NA"
			tmp1 = outs[mypep][myid].split(";")
			pfam_info = pfam_info + outs[mypep][myid] + ";"
			for item in tmp1:
				tmp2 = item.split(":")
				ann1 = "NA"
				ann2 = "NA"
				if tmp2[0] in pfams:
					ann1 = pfams[tmp2[0]]
				if tmp2[1] in pfams:
					ann2 = pfams[tmp2[1]]
				ann3 = ann1 + ":" + ann2
				if myann == "NA":
					myann = ann3
				else:
					myann = myann + ";" + ann3
			ann_info = ann_info + myann + ";"
			level_info = level_info + mylevel + ";"
		# foreach item
		pfam_info = re.sub(";$", "", pfam_info)
		ann_info = re.sub(";$", "", ann_info)
		level_info = re.sub(";$", "", level_info)
		if pfam_info == "":
			continue
		open_out.write(mypep + "\tDOMINE_interaction" + "\t" + pfam_info + "\t" + ann_info + "\t" + level_info + "\n")
	# foreach peptide family
	open_out.close()
# assign_annotation


#==============================================================
# annotate DDIs
#==============================================================
def DDI_annotation (extension, pfam_path, filter_flag, spe_level, interact, pfams, human_pfam, suffix):
	filelist = utilities.find_files(pfam_path, extension, None)
	for myfile in filelist:
		#myfile = pfam_path + "/" + samplelist + "/" + samplelist + ".interpro.PfamDomain.tsv"
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			#print ("OK!\t" + myfile)
			myout = re.sub(extension, suffix, myfile)
			myout_detail = re.sub(".tsv", ".detail.tsv", myout)
			peptide = collect_pfam_info (myfile)
			assign_interaction (filter_flag, spe_level, pfams, human_pfam, interact, peptide, myout, myout_detail)
        # if exist file
    # foreach split part
# DDI_annotation


#==============================================================
###########  Main processing ############
#=============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start ddi_DOMINE_protein.py -p " + values.path + " ####\n")


	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	pfams = collect_pfam_ann (config.pfam_database)
	human_pfam = collect_mapping_info (config.human_pfam_database)
	interact = collect_interaction_info (config.interaction_database, values.filter, human_pfam)
	sys.stderr.write("Get info ......done\n")
	
	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign interaction to peptide families ......starting\n")
	DDI_annotation (values.extension, values.path, values.filter, "HC", interact, pfams, human_pfam, values.suffix)
	sys.stderr.write("\nAssign interaction to peptide families ......done\n")

	sys.stderr.write("### Finish ddi_DOMINE_protein.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

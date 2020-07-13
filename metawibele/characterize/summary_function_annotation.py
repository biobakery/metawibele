#!/usr/bin/env python

"""
MetaWIBELE: summary_function_annotation module
Summary function annotations

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
Summary functional annotation
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-l', "--list",
	                    help='input the list file of annotation files',
	                    required=True)
	parser.add_argument('-a', "--uniref-annotation",
	                    help='input uniref90 annotation file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output functional annotation file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_cluster_info (clust_file):	# discovery_cohort.peptides.clust
	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	myclust_id = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			mym = re.search("cluster=([\d]+)", line)
			myclust_id = "Cluster_" + mym.group(1)
			if not myclust in cluster:
				cluster[myclust] = {}
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myclust][myid] = myclust_id
	# foreach line
	open_file.close()
	return cluster
# function collect_cluster_info


#==============================================================
# collect UniRef annotation info
#==============================================================
def collect_uniref_info (clust_file):   # summary_peptide_family_annotation.uniref90_annotation.tsv
	note = {}
	titles = {}
	open_file = open(clust_file, "r")
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
		mynote = info[titles["note"]]
		note[info[0]] = mynote
	# foreach line
	open_file.close()
	return note, id_flag
# function collect_uniref_info


#==============================================================
# assign annotation to peptide families info
#==============================================================
def collect_annotation (list_file, id_flag):
	annotation = {}
	anns = {}
	note = {}
	open_list = open(list_file, "r")
	for myfile in open_list.readlines():
		myfile = myfile.strip()
		if not len(myfile):
			continue
		if re.search("^#", myfile):
			continue
		if not os.path.isfile(myfile):
			print("File not exist!\t" + myfile)
			continue
		open_file = open(myfile, "r")
		titles = {}
		mym = re.search(config.basename + "_([\S]+)_proteinfamilies", os.path.basename(myfile))
		method = mym.group(1)
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			if re.search("^" + id_flag, line):
				for item in info:
					titles[item] = info.index(item)
				continue
			myid = info[titles[id_flag]]
			if "note" in titles:
				note[myid] = info[titles["note"]]
			if re.search("UniRef", method):
				desc = info[titles["Protein_names"]]
				tax = info[titles["Tax"]]
				taxID = info[titles["TaxID"]]
				reptax = info[titles["Rep_Tax"]]
				reptaxID = info[titles["Rep_TaxID"]]
				#org = info[titles["organism"]]
				uniprot = info[titles["UniProtKB"]]
				uniref = info[titles["unirefID"]]
				if len(info) < 8:
					print(line)
				mytype = info[1] + "\t" + info[2] + "\t" + desc + "\t" + tax + "\t" + taxID + "\t" +  reptax + "\t" + reptaxID + "\t" + uniprot + "\t" + uniref
			else:
				mytype = info[1] + "\t" + info[2] + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
			if not myid in annotation:
				annotation[myid] = {}
			if not method in annotation[myid]:
				annotation[myid][method] = []
			annotation[myid][method].append(mytype)
			if not myid in anns:
				anns[myid] = {}
			if not method in anns[myid]:
				anns[myid][method] = {}
			anns[myid][method][info[1]] = ""
		# foreach line
		open_file.close()
	# foreach annotation file
	open_list.close()

	# collect clusters which have on decent homologies in UniRef90 or uncharacterized in UniRef90
	annotation_non_uniref = {}
	anns_uniref = {}
	for myid in annotation.keys():
		if "UniRef90" in anns[myid]:
			flag = 0
			if "UniRef90_unknown" in anns[myid]["UniRef90"] or "UniRef90_uncharacterized" in anns[myid]["UniRef90"]:	# UniRef90 unannotated ones
				#print(myid)
				flag = 1
				if not myid in annotation_non_uniref:
					annotation_non_uniref[myid] = {}
				for method in annotation[myid]:
					if not method in annotation_non_uniref[myid]:
						annotation_non_uniref[myid][method] = []
					for item in annotation[myid][method]:
						annotation_non_uniref[myid][method].append(item)
				# foreach method
			# if unannotated in UniRef90
			if flag == 0: # characterized by UniRef90	
				if not myid in anns_uniref:
					anns_uniref[myid] = {}
				for method in annotation[myid]:
					if not method in anns_uniref[myid]:
						anns_uniref[myid][method] = []
					for item in annotation[myid][method]:
						anns_uniref[myid][method].append(item)
				# foreach method
		# if have UniRef method
	# foreach cluster
	return annotation, anns_uniref, annotation_non_uniref, note
# collect_annotation

# assign_annotation
def assign_annotation (id_flag, pep_cluster, annotation, study, note1, note2, outfile):
	outs = {}
	anns = {}
	number = {}
	for pepid in sorted(pep_cluster.keys()): # foreach peptide family
		flag = 0
		for member in pep_cluster[pepid].keys():
			if id_flag == utilities.PROTEIN_FAMILY_ID:
				if pepid != member:
					continue
				myclust_id = pep_cluster[pepid][member]
			if id_flag == utilities.PROTEIN_ID:
				myclust_id = member
			if not myclust_id in annotation:	# no corresponding annotation
				continue
			if not myclust_id in outs:
				outs[myclust_id] = {}
			for method in sorted(annotation[myclust_id].keys()):
				if not method in outs[myclust_id]:
					outs[myclust_id][method] = {}
				for myid in annotation[myclust_id][method]:
					flag = 1
					myinfo = myid.split("\t")
					mytype = myinfo[0]
					outs[myclust_id][method][myid] = ""
					if not myclust_id in anns:
						anns[myclust_id] = {}
					anns[myclust_id][mytype] = ""
			# foreach type
		# foreach member
	# foreach peptide cluster

	interest = ["cellWall", "outerMembrane", "extracellular", "signaling", "transmembrane", "interaction", "PfamDomain", "COG", "KEGG-KOs", "GO"]
	other = ["cytoplasmic", "cytoplasmicMembrane", "periplasmic",
	         "cellSurface", "cytoplasm", "membrane", "periplasm", "nucleus", "fimbrium", "virion", "sporeCore", "mitochondrion", "cellEnvelop", "cellMembrane", "cellInnerMembrane", "bacterialFlagellum",
	         "others", "Others"]
	unknown = ["unknown", "hypothetical"]
	uncharacterized = ["uncharacterized"]
	interpro = ["SUPERFAMILY", "ProSitePatterns", "ProSiteProfiles", "Gene3D", "PANTHER", "TIGRFAM", "SFLD", "ProDom", "Hamap", "SMART", "CDD", "PRINTS", "PIRSF", "MobiDBLite", "Coils"]
	total_num = 0
	for myclust in anns.keys():
		total_num = total_num + 1
		interest_flag = 0
		other_flag = 0
		unchar_flag = 0
		pfam_flag = 0
		interpro_flag = 0
		for mytype in anns[myclust].keys():
			if not re.search("_", mytype):
				# debug
				print(mytype)
				continue
			mym, category = mytype.split("_")
			category = re.sub("GO\(BP\)", "GO", category)
			category = re.sub("GO\(CC\)", "GO", category)
			category = re.sub("GO\(MF\)", "GO", category)
			if category in interest:
				interest_flag = 1
				if not category in number:
					number[category] = {}
				number[category][myclust] = ""
			if category in other or category in interpro:
				other_flag = 1
			if category in uncharacterized:
				unchar_flag = 1
			#if category in interpro:
			#	interpro_flag = 1
		# foreach type of annotation
		if interest_flag == 0 and other_flag == 1:	# other type annotation
			category = "Others"	
			if not category in number:
				number[category] = {}
			number[category][myclust] = ""
		if interest_flag == 0 and other_flag == 0 and unchar_flag == 1:	# uncharacterized type annotation
			category = "uncharacterized"	
			if not category in number:
				number[category] = {}
			number[category][myclust] = ""
		if interest_flag == 0 and other_flag == 0 and unchar_flag == 0:	# unknown type annotation
			category = "Unknown"	
			if not category in number:
				number[category] = {}
			number[category][myclust] = ""
	# foreach cluster

	open_out = open(outfile, "w")
	open_out.write(id_flag + "\tstudy\tmethod\tcategory\ttype\tdetail\tProtein_names\tTax\tTaxID\tRep_Tax\tRep_TaxID\tUniProtKB\tunirefID\tnote\n")
	for myclust in sorted(outs.keys()):
		note = {}
		if myclust in note1:
			note[note1[myclust]] = ""
		if myclust in note2:
			note[note2[myclust]] = ""
		mynote = "good"
		if len(note.keys()) > 0:
			mynote = ";".join(sorted(note.keys()))
			#if len(note.keys()) > 1:
			#	if "good" in note:
			#		note.pop("good") 
			#	mynote = ";".join(sorted(note.keys()))
		for method in sorted(outs[myclust].keys()):
			for myid in sorted(outs[myclust][method].keys()):
				myinfo = myid.split("\t")
				mytype = myinfo[0]
				if not re.search("_", mytype):
					# debug
					print(mytype)
					continue
				mym, category = mytype.split("_")
				if category in unknown:
					category = "Unknown"
				if category in other:
					category = "Others"
				open_out.write(myclust + "\t" + study + "\t" + method + "\t" + category + "\t" + myid + "\t" + mynote + "\n")
			# foreach type
		# foreach method
	# foreach peptide family
	open_out.close()
# assign_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start summary_function_annotation.py -l " + values.list + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	pep_cluster = collect_cluster_info (config.protein_family)
	sys.stderr.write("Get cluster info ......done\n")
	
	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	note1, id_flag = collect_uniref_info (values.uniref_annotation)
	annotation, anns_uniref, anns_non_uniref, note2 = collect_annotation (values.list, id_flag)
	sys.stderr.write("Get annotation info ......done\n")

	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign annotation......starting\n")
	assign_annotation (id_flag, pep_cluster, annotation, config.study, note1, note2, values.output)
	uniref_out = re.sub(".tsv", ".uniref.tsv", values.output)
	assign_annotation (id_flag, pep_cluster, anns_uniref, config.study, note1, note2, uniref_out)
	uniref_non = re.sub(".tsv", ".non_uniref.tsv", values.output)
	assign_annotation (id_flag, pep_cluster, anns_non_uniref, config.study, note1, note2, uniref_non)
	sys.stderr.write("\nAssign annotation......done\n")

	sys.stderr.write("### Finish summary_function_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

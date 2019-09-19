#!/usr/bin/env python

"""
MetaWIBELE: summary_protein_annotation module
Combine annotation for each ORF

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
Combine annotation for each ORF
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-l', help='input the list file of annotation file', required=True)
	parser.add_argument('-o', help='output cluster distribution file', required=True)
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
# assign annotation to peptide families info
#==============================================================
def collect_annotation(list_file):
	annotation = {}
	anns = {}
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
		mym = re.search("summary_([\S]+)_protein_family", myfile)
		method = mym.group(1)
		titles = {}
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			if re.search("^seqID", line):
				for item in info:
					titles[item] = info.index(item)
				continue
			myid = info[0]
			mytype1 = "NA"
			mytype2 = "NA"
			mytype3 = "NA"
			mytype4 = "NA"
			mytype5 = "NA"
			if re.search("UniRef", method):
				go1 = info[titles["GO(BP)"]]
				go2 = info[titles["GO(MF)"]]
				go3 = info[titles["GO(CC)"]]
				kegg = info[titles["KO"]]
				cog = info[titles["COG"]]
				tax = info[titles["Tax"]]
				taxID = info[titles["TaxID"]]
				reptax = info[titles["Rep_Tax"]]
				reptaxID = info[titles["Rep_TaxID"]]
				org = info[titles["organism"]]
				uniprot = info[titles["UniProtKB"]]
				uniref = info[titles["unirefID"]]	
				if len(info) < 8:
					print(line)
				#mytype = info[1] + "\t" + info[2] + "\t" + info[3] + "\t" + info[4] + "\t" + info[5] + "\t" + info[6] + "\t" + info[7]
				mytype = info[1] + "\t" + info[2] + "\t" + tax + "\t" + taxID + "\t" +  reptax + "\t" + reptaxID + "\t" + org + "\t" + uniprot + "\t" + uniref
				if go1 != "NA":
					mytype1 = "UniRef90_GO(BP)" + "\t" + go1 + "\t" + tax + "\t" + taxID + "\t" + reptax + "\t" + reptaxID + "\t" + org + "\t" + uniprot + "\t" + uniref 
				if go2 != "NA":
					mytype2 = "UniRef90_GO(MF)" + "\t" + go2 + "\t" + tax + "\t" + taxID + "\t" + reptax + "\t" + reptaxID + "\t" + org + "\t" + uniprot + "\t" + uniref
				if go3 != "NA":
					mytype3 = "UniRef90_GO(CC)" + "\t" + go3 + "\t" + tax + "\t" + taxID + "\t" + reptax + "\t" + reptaxID + "\t" + org + "\t" + uniprot + "\t" + uniref 
				if kegg != "NA":
					mytype4 = "UniRef90_KEGG-KOs" + "\t" + kegg + "\t" + tax + "\t" + taxID + "\t" + reptax + "\t" + reptaxID + "\t" + org + "\t" + uniprot + "\t" + uniref
				if cog != "NA":
					mytype5 = "UniRef90_COG" + "\t" + cog + "\t" + tax + "\t" + taxID + "\t" + reptax + "\t" + reptaxID + "\t" + org + "\t" + uniprot + "\t" + uniref
			else:
				mytype = info[1] + "\t" + info[2] + "\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
			if not myid in annotation:
				annotation[myid] = {}
			if not method in annotation[myid]:
				annotation[myid][method] = []
			annotation[myid][method].append(mytype)
			if mytype1 != "NA":
				annotation[myid][method].append(mytype1)
			if mytype2 != "NA":
				annotation[myid][method].append(mytype2)
			if mytype3 != "NA":
				annotation[myid][method].append(mytype3)
			if mytype4 != "NA":
				annotation[myid][method].append(mytype4)
			if mytype5 != "NA":
				annotation[myid][method].append(mytype5)
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
	return annotation, anns_uniref, annotation_non_uniref
# collect_annotation

# assign_annotation
def assign_annotation (pep_cluster, annotation, study, outfile):
	outs = {}
	anns = {}
	number = {}
	for pepid in sorted(pep_cluster.keys()): # foreach peptide family
		flag = 0
		for member in pep_cluster[pepid].keys():
			#if pepid != member:
			#	continue
			#myclust_id = pep_cluster[pepid][member]
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

	#interest = ["secreted", "cellWall", "outerMembrane", "extracellular", "signaling", "transmembrane", "uncharacterized", "interaction", "cellSurface"]
	interest = ["secreted", "cellWall", "outerMembrane", "extracellular", "signaling", "transmembrane", "interaction", "PfamDomain", "SUPERFAMILY"]
	other = ["cytoplasmic", "cytoplasmicMembrane", "periplasmic", "others", "Others", "GO(BP)", "GO(CC)", "GO(MF)", "KEGG-KOs", "COG"]
	unknown = ["unknown", "hypothetical"]
	uncharacterized = ["uncharacterized"]
	interpro = ["ProSitePatterns", "ProSiteProfiles", "Gene3D", "PANTHER", "TIGRFAM", "SFLD", "ProDom", "Hamap", "SMART", "CDD", "PRINTS", "PIRSF", "MobiDBLite", "Coils"]
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
				#print(mytype)
				continue
			mym, category = mytype.split("_")
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
	open_out.write("seqID\tstudy\tmethod\tcategory\ttype\tdetail\tTax\tTaxID\tRep_Tax\tRep_TaxID\torganism\tUniProtKB\tunirefID\n")
	for myclust in sorted(outs.keys()):
		for method in sorted(outs[myclust].keys()):
			for myid in sorted(outs[myclust][method].keys()):
				myinfo = myid.split("\t")
				mytype = myinfo[0]
				if not re.search("_", mytype):
					# debug
					#print(mytype)
					continue
				mym, category = mytype.split("_")
				if category in unknown:
					category = "Unknown"
				if category in other:
					category = "Others"
				open_out.write(myclust + "\t" + study + "\t" + method + "\t" + category + "\t" + myid + "\n")
			# foreach type
		# foreach method
	# foreach peptide family
	open_out.close()
	
#	outfile2 = re.sub(".tsv", ".plot.tsv", outfile)
#	open_out = open(outfile2, "w")
#	open_out.write("category\tnumber\tpercentage\n")
#	for category in sorted(number.keys()):
#		mynum = len(number[category].keys())
#		myper = float(mynum)/float(total_num)
#		open_out.write(category + "\t" + str(mynum) + "\t" + str(myper) + "\n")
#	open_out.close()
# assign_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start summary_protein_annotation.py -l " + values.l + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	pep_cluster = collect_cluster_info (config.protein_family)
	sys.stderr.write("Get cluster info ......done\n")
	
	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	annotation, anns_uniref, anns_non_uniref = collect_annotation(values.l)
	sys.stderr.write("Get annotation info ......done\n")

	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign annotation to peptide families ......starting\n")
	assign_annotation (pep_cluster, annotation, config.study, values.o)
	uniref_out = re.sub(".tsv", ".uniref.tsv", values.o)
	assign_annotation (pep_cluster, anns_uniref, config.study, uniref_out)
	uniref_non = re.sub(".tsv", ".non_uniref.tsv", values.o)
	assign_annotation (pep_cluster, anns_non_uniref, config.study, uniref_non)
	sys.stderr.write("\nAssign annotation to peptide families ......done\n")

	sys.stderr.write("### Finish summary_protein_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

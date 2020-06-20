#!/usr/bin/env python

"""
MetaWIBELE: uniref_protein_family module
Summary the uniref annotation results by mapping protein clusters to UniRef

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
Summary the uniref annotation results from mapping protein clusters to UniRef
"""

def get_args ():
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument("-u", "--protein-annotation",
	                   help='input uniref annotated file',
	                   required=True)
	parser.add_argument("-m", "--mapping",
	                    help='input uniref mapping file',
	                    required=True)
	parser.add_argument("-f", "--method",
	                    help='specify how to assign annotations to families',
	                    choices=["centroid", "consistency"],
	                    required=True,
	                    default="centroid")
	parser.add_argument("-o", "--output",
	                    help='output annotation detailed file for UniRef annotation',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_peptide_cluster_info (clust_file):	# discovery_cohort.peptides.clust
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
# function collect_peptide_cluster_info

def collect_gene_cluster_info (clust_file):	# discovery_cohort.genes.clust
	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myclust
	# foreach line
	open_file.close()
	return cluster
# function collect_gene_cluster_info


#==============================================================
# collect uniref mapping info
#==============================================================
def collect_uniref_mapping (mapfile):
	maps = {}
	titles = {}
	open_file = open(mapfile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search(utilities.PROTEIN_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myid = info[titles[utilities.PROTEIN_ID]]
		mytype1 = info[titles["query_type"]]
		mytype2 = info[titles["mutual_type"]]
		myiden = info[titles["identity"]]
		mycov1 = info[titles["query_coverage"]]
		mycov2 = info[titles["mutual_coverage"]]
		maps[myid] = mytype1 + "\t" + mytype2 + "\t" + myiden + "\t" + mycov1 + "\t" + mycov2
	# foreach line
	open_file.close()
	return maps
# collect_uniref_mapping


#==============================================================
# assign annotation to peptide families info
#==============================================================
def collect_annotation (uniref): 
	locations = ["Bacterial flagellum",
				 "cell wall",
	             "Cell envelop", "Cell membrane", "Cell inner membrane", "Cell outer membrane", "Cell surface",
	             "Fimbrium",
	             "Nucleus",
	             "Periplasm",
	             "Mitochondrion", "Mitochondrion inner membrane",
	             "Secreted",
	             "Spore core",
	             "Virion",
	             "Membrane",
	             "Cytoplasm"]
	annotation = {}
	titles = {}
	open_file = open(uniref, "r")
	line = open_file.readline()
	line = line.strip()
	ann_title = line
	ann_title = re.sub("^name\t", "", ann_title) 
	info = line.split("\t")
	myindex = 0
	while myindex < len(info):
		item = info[myindex]
		titles[item] = myindex
		myindex = myindex + 1
	unitype = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if myid in annotation:
			continue
		mytype = info[1]
		myinfo = "\t".join(info[1:len(info)])
		ann = info[titles["Protein_names"]]
		if ann == "NA":
			ann = mytype
		if not re.search("_unknown", info[1]):
			mym = re.search("^([^_]+)", info[1])
			unitype = mym.group(1)
			#go_bp = info[titles["GO(BP)"]]
			#go_mf = info[titles["GO(MF)"]]
			#go_cc = info[titles["GO(CC)"]]
			if "GO" in titles:
				go = info[titles["GO"]]
			else:
				go = "NA"
			if "GO_BP" in titles:
				go_bp = info[titles["GO_BP"]]
			else:
				go_bp = "NA"
			if "GO_cc" in titles:
				go_cc = info[titles["GO_CC"]]
			else:
				go_cc = "NA"
			if "GO_MF" in titles:
				go_MF = info[titles["GO_MF"]]
			else:
				go_mf = "NA"
			if "KO" in titles:
				kegg = info[titles["KO"]]
			else:
				kegg = "NA"
			if "eggNOG" in titles:
				cog = info[titles["eggNOG"]]
			else:
				cog = "NA"
			if "Pfam" in titles:
				pfam = info[titles["Pfam"]]
			else:
				pfam = "NA"
			if "Level4EC" in titles:
				level4ec = info[titles["Level4EC"]]
			else:
				level4ec = "NA"
			interest_flag = 0

			# based on GO info to define uncharacterized protein
			if go == "NA" and go_bp == "NA" and go_cc == "NA" and go_mf == "NA":
				mytype = unitype + "_uncharacterized"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + mytype + "\t" + ann + "\n" + myinfo)

			# functional annotation
			if go != "NA":
				interest_flag = 1
				mytype = unitype + "_GO"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + go + "\t" + ann + "\n" + myinfo)
			if go_bp != "NA":
				interest_flag = 1
				mytype = unitype + "_GO-BP"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + go_bp + "\t" + ann + "\n" + myinfo)
			if go_cc != "NA":
				interest_flag = 1
				mytype = unitype + "_GO-CC"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + go_cc + "\t" + ann + "\n" + myinfo)
			if go_mf != "NA":
				interest_flag = 1
				mytype = unitype + "_GO-MF"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + go_mf + "\t" + ann + "\n" + myinfo)
			if kegg != "NA":
				interest_flag = 1
				mytype = unitype + "_KEGG-KOs"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + kegg + "\t" + ann + "\n" + myinfo)
			if cog != "NA":
				interest_flag = 1
				mytype = unitype + "_eggNOG"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + cog + "\t" + ann + "\n" + myinfo)
			if pfam != "NA":
				interest_flag = 1
				mytype = unitype + "_PfamDomain"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + pfam + "\t" + ann + "\n" + myinfo)
			if level4ec != "NA":
				interest_flag = 1
				mytype = unitype + "_Level4EC"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + level4ec + "\t" + ann + "\n" + myinfo)
		# if not unknown
		else:
			if not myid in annotation:
				annotation[myid] = []
			annotation[myid].append(mytype + "\t" + mytype + "\t" + ann + "\n" + myinfo)
	# foreach line
	open_file.close()
	return annotation, unitype, ann_title
# collect_annotation

# assign_annotation
def assign_annotation (identity_cutoff, coverage_cutoff, cutoff, pep_cluster, annotation, unitype, ann_title, maps, assign_flag, outfile_detail):
	outs = {}
	outs_ORF = {}
	outs_info = {}
	number = {}
	anns = {}
	titles = {}
	maps_num = {}
	pers = {}
	ann_info = {}
	outs2 = {}
	outs2_info = {}
	consistency = {}
	cluster = {}
	tmp = ann_title.split("\t")
	for item in tmp:
		titles[item] = tmp.index(item)
	for pepid in sorted(pep_cluster.keys()): # foreach peptide family
		flag = 0
		detail_info = {}
		mytotal = 0
		clust_id = ""
		if pepid in maps:
			myquery, mymutual, myiden, query_cov, mutual_cov = maps[pepid].split("\t")
			myt = mymutual
			#if mymutual == "low_confidence" and myquery == "high_confidence":
			#	myt = "low_confidence_high_query"
			#if mymutual == "low_confidence" and myquery == "low_confidence":
			#	myt = "low_confidence_others"
			if myt == "high_confidence":
				myt = "UniRef90_strong_homology"
			if myt == "low_confidence":
				if float(myiden) < float(identity_cutoff) or float(mutual_cov) < float(coverage_cutoff):
					myt = "UniRef90_worse_homology"
				else:
					myt = "UniRef90_weak_homology"
			if myt == "no_hit":
				myt = "UniRef90_worse_homology"
			if not myt in maps_num:
				maps_num[myt] = {}
			maps_num[myt][pepid] = ""
		else:
			# debug
			print("No UniRef90 mapping info!\t" + pepid)
			myt = "UniRef90_worse_homology"
			myp = "NA"
			if not myt in maps_num:
				maps_num[myt] = {}
			maps_num[myt][pepid] = ""
		for member in pep_cluster[pepid].keys():
			#if member != pepid:	# use the annotation of representative
			#	continue
			#if not member in gene_cluster:	# no corresponding gene cluster
			#	print("Peptide ID has no corresponding gene cluster!\t" + member)
			#	continue
			#gene_id = gene_cluster[member]
			gene_id = member
			clust_id = pep_cluster[pepid][member]
			if not gene_id in annotation:	# no corresponding annotation
				# debug
				print("No UniRef90 mapping info for the member!\t" + gene_id)
				continue
			for myall in annotation[gene_id]:
				myid, myinfo = myall.split("\n")
				mytype, mytype1, ann = myid.split("\t")
				tmp_info = myinfo.split("\t")
				if re.search("_unknown", tmp_info[0]):
					tmp_info[0] = "NA"
				mygene_name = "NA"
				mytax = "NA"
				mytaxID = "NA"
				myreptax = "NA"
				myreptaxID = "NA"
				myorg = "NA"
				myuniprot = "NA"
				myunref = "NA"
				if "Gene_names" in titles:
					mygene_name = tmp_info[titles["Gene_names"]]
				if "Tax" in titles:
					mytax = tmp_info[titles["Tax"]]
				if "TaxID" in titles:
					mytaxID = tmp_info[titles["TaxID"]]
				if "Rep_Tax" in titles:
					myreptax = tmp_info[titles["Rep_Tax"]]
				if "Rep_TaxID" in titles:
					myreptaxID = tmp_info[titles["Rep_TaxID"]]
				#if "Organism" in titles:
				#	myorg = tmp_info[titles["Organism"]]
				if "UniProtKB" in titles:
					myuniprot= tmp_info[titles["UniProtKB"]]
				if "UniRefID" in titles:
					myuniref = tmp_info[titles["UniRefID"]]
				if member == pepid:	# use the annotation of representative
					flag = 1
					if not clust_id in outs:
						outs[clust_id] = {}
					if not clust_id in outs_info:
						outs_info[clust_id] = {}
					outs[clust_id][mytype + "\t" + mytype1 + "\t" + ann + "\t" + mygene_name + "\t" + mytax + "\t" + mytaxID + "\t" + myreptax + "\t" + myreptaxID + "\t" + myuniprot + "\t" + myuniref] = ""
					outs_info[clust_id][mytype + "\t" + myinfo] = ""
					if not clust_id in anns:
						anns[clust_id] = {}
					anns[clust_id][mytype] = ""
				# for cluster
				if not member in outs_ORF:
					outs_ORF[member] = {}
				#outs_ORF[member][mytype + "\t" + mytype1 + "\t" + ann + "\t" + tmp_info[titles["Gene_names"]] + "\t" + tmp_info[titles["Tax"]] + "\t" + tmp_info[titles["TaxID"]] + "\t" + tmp_info[titles["Rep_Tax"]] + "\t" + tmp_info[titles["Rep_TaxID"]] + "\t" + tmp_info[titles["Organism"]] + "\t" + tmp_info[titles["UniProtKB"]] + "\t" + tmp_info[titles["UniRefID"]]] = ""
				outs_ORF[member][mytype + "\t" + mytype1 + "\t" + ann + "\t" + mygene_name + "\t" + mytax + "\t" + mytaxID + "\t" + myreptax + "\t" + myreptaxID  + "\t" + myuniprot + "\t" + myuniref] = ""
				if tmp_info[0] != "NA":
					mytotal = mytotal + 1
				myuniref_id = tmp_info[0]
				if not clust_id in outs2_info:
					outs2_info[clust_id] = {}
				if not member in outs2_info[clust_id]:
					outs2_info[clust_id][member] = myuniref_id
				if not myuniref_id in detail_info:
					detail_info[myuniref_id] = {}
				detail_info[myuniref_id][member] = ""
				if not myuniref_id in ann_info:
					ann_info[myuniref_id] = {}
				#ann_info[myuniref_id][mytype + "\t" + mytype1 + "\t" + ann + "\t" + tmp_info[titles["Gene_names"]] + "\t" + tmp_info[titles["Tax"]] + "\t" + tmp_info[titles["TaxID"]] + "\t" + tmp_info[titles["Rep_Tax"]] + "\t" + tmp_info[titles["Rep_TaxID"]] + "\t" + tmp_info[titles["Organism"]] + "\t" + tmp_info[titles["UniProtKB"]] + "\t" + tmp_info[titles["UniRefID"]]] = ""
				ann_info[myuniref_id][mytype + "\t" + mytype1 + "\t" + ann + "\t" + mygene_name + "\t" + mytax + "\t" + mytaxID + "\t" + myreptax + "\t" + myreptaxID + "\t" + myuniprot + "\t" + myuniref] = ""
			# foreach type
		# foreach member
		if flag == 0:
			for member in pep_cluster[pepid].keys():
				#if member != pepid:	# use the annotation of representative
				#	continue
				clust_id = pep_cluster[pepid][member]
			if not clust_id in outs:
				outs[clust_id] = {}
				outs_info[clust_id] = {}
				type_tmp = unitype + "_unknown"
				mystr = type_tmp
				tmp1 = ann_title.split("\t")
				item_num = len(tmp1)
				mynum = 1
				while mynum <= item_num:
					mystr = mystr + "\tNA"
					mynum = mynum + 1
				mystr2 = "NA\tNA\tNA\tNA\tNA\tNA\tNA"
				outs[clust_id][type_tmp + "\t" + type_tmp + "\t" + type_tmp  + "\t" + mystr2] = ""
				outs_info[clust_id][mystr] = ""
				if not clust_id in anns:
					anns[clust_id] = {}
				anns[clust_id][type_tmp] = ""
				myuniref_id = "NA"
				for member in pep_cluster[pepid].keys():
					if not member in outs_ORF:
						outs_ORF[member] = {}
					outs_ORF[member][type_tmp + "\t" + type_tmp + "\t" + type_tmp + "\t" + mystr2] = ""
					if not myuniref_id in detail_info:
						detail_info[myuniref_id] = {}
					detail_info[myuniref_id][member] = ""
				if not myuniref_id in ann_info:
					ann_info[myuniref_id] = {}
				ann_info[myuniref_id][type_tmp + "\t" + type_tmp + "\t" + type_tmp + "\t" + mystr2] = ""
		# if no annotation
		cluster[clust_id] = ""

		# consistency
		detail = {}
		mytotal = len(pep_cluster[pepid].keys())
		for mytype in detail_info.keys():
			mynum = len(detail_info[mytype].keys())
			if not mynum in detail:
				detail[mynum] = []
			detail[mynum].append(mytype)
		# foreach type
		for mynum in sorted(detail.keys(), key=int, reverse=True):
			item_num = len(detail[mynum])
			for item in detail[mynum]:
				myper = round(float(mynum)/float(mytotal), 2)
				if mynum != 0:
					if not myper in pers:
						pers[myper] = {}
					if not mytotal in pers[myper]:
						pers[myper][mytotal] = {}
					pers[myper][mytotal][clust_id] = ""
					if myper >= float(cutoff):
						if not clust_id in outs2:
							outs2[clust_id] = {}
						if item in ann_info:
							for myann in ann_info[item]:
								outs2[clust_id][myann] = ""
						mylabel = item
						myrep = "NA"
						if clust_id in outs:
							for myann in outs[clust_id].keys():
								tmp = myann.split("\t")
								myrep = tmp[8]	# the UniRefID
						myagree = "no"
						if mylabel == myrep:
							myagree = "yes"
						consistency[clust_id] = mylabel + "\t" + myrep + "\t" + myagree
					else:	# less consistency using representative results
						if clust_id in outs:
							if not clust_id in outs2:
								outs2[clust_id] = {}
							for myann in outs[clust_id].keys():
								outs2[clust_id][myann] = ""
			# foreach type
			break
		# foreach number
	# foreach peptide cluster

	#### output info ####
	outfile1 = re.sub(".detail.tsv", ".ORF.detail.tsv", outfile_detail)
	open_out = open(outfile_detail, "w")
	open_out1 = open(outfile1, "w")
	tmp1 = ann_title.split("\t")
	open_out.write(utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tProtein_names\tGene_names\tTax\tTaxID\tRep_Tax\tRep_TaxID\tUniProtKB\tunirefID" + "\n")
	open_out1.write(utilities.PROTEIN_ID + "\ttype\tdetail\tProtein_names\tGene_names\tTax\tTaxID\tRep_Tax\tRep_TaxID\tUniProtKB\tunirefID" + "\n")
	for myclust in sorted(cluster.keys()):
		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in outs:
				for myid in sorted(outs[myclust].keys()):
					open_out.write(myclust + "\t" + myid + "\n")
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in outs2:
				for myid in sorted(outs2[myclust].keys()):
					open_out.write(myclust + "\t" + myid + "\n")
	# foreach cluster
	open_out.close()
	for mypep in sorted(outs_ORF.keys()):
		for myid in sorted(outs_ORF[mypep].keys()):
			open_out1.write(mypep + "\t" + myid + "\n")
		# foreach type
	# foreach peptide family
	open_out1.close()
	
# assign_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start uniref_protein_family.py -m " + values.mapping + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	pep_cluster = collect_peptide_cluster_info (config.protein_family)
	sys.stderr.write("Get cluster info ......done\n")
	
	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	maps = collect_uniref_mapping(values.mapping)
	annotation, uniref_type, ann_title = collect_annotation(values.protein_annotation)
	sys.stderr.write("Get annotation info ......done")

	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign annotation to peptide families ......starting\n")
	assign_annotation (config.tshld_identity, config.tshld_coverage, config.tshld_consistency, pep_cluster, annotation, uniref_type, ann_title, maps, values.method, values.output)
	sys.stderr.write("\nAssign annotation to peptide families ......done\n")

	sys.stderr.write("### Finish uniref_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python3

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
def collect_uniref_mapping (mapfile):	# discovery_cohort_meta.genes.clust.rep.uniref90.stat.tsv
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
		if re.search("^name", line):
			continue
		info = line.split("\t")
		myid = info[0]
		if myid in annotation:
			continue
		mytype = info[1]
		myinfo = "\t".join(info[1:len(info)])
		ann = info[titles["Description"]]
		if ann == "NA":
			ann = mytype
		if not re.search("_unknown", info[1]):
			mym = re.search("^([^_]+)", info[1])
			unitype = mym.group(1)
			go_bp = info[titles["GO(BP)"]]
			go_mf = info[titles["GO(MF)"]]
			go_cc = info[titles["GO(CC)"]]
			kegg = info[titles["KO"]]
			cog = info[titles["COG"]]
			pfam = info[titles["Pfam"]]
			transmembrane = info[titles["Transmembrane"]]
			signaling = info[titles["Signal_peptide"]]
			secreted = info[titles["Subcellular_location"]]
			if re.search("^NA ", secreted):
				secreted = "NA"
			interest_flag = 0

			# based on GO info to define uncharacterized protein
			#if info[titles["GO(BP)"]] == "NA" and info[titles["GO(MF)"]] == "NA" and info[titles["GO(CC)"]] == "NA" and info[titles["Pfam"]] == "NA":
			if info[titles["GO(BP)"]] == "NA" and info[titles["GO(MF)"]] == "NA" and info[titles["GO(CC)"]] == "NA":
				mytype = unitype + "_uncharacterized"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + mytype + "\t" + ann + "\n" + myinfo)

			# secondary structure
			if transmembrane != "NA":
				interest_flag = 1
				tmp = transmembrane.split(";")
				tmp = list(dict.fromkeys(tmp))
				transmembrane = ";".join(sorted(tmp))
				mytype = unitype + "_transmembrane"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + transmembrane + "\t" + ann + "\n" + myinfo)
			if signaling != "NA":
				interest_flag = 1
				tmp = signaling.split(";")
				tmp = list(dict.fromkeys(tmp))
				signaling = ";".join(sorted(tmp))
				mytype = unitype + "_signaling"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + signaling + "\t" + ann + "\n" + myinfo)
			if secreted != "NA":
				myflag = 0
				if re.search("extracellular", go_cc) or re.search("Extracellular", go_cc):
					if re.search("extracellular", secreted) or re.search("Extracellular", secreted):
						myflag = 1
						interest_flag = 1
						mytype = unitype + "_extracellular"
						if not myid in annotation:
							annotation[myid] = []
						annotation[myid].append(mytype + "\t" + secreted + "\t" + ann + "\n" + myinfo)
					else:
						if re.search("secretion", secreted) or re.search("secretory", secreted) or re.search("secreted", secreted) or re.search("secretin", secreted) or re.search("Secretion", secreted) or re.search("Secretory", secreted) or re.search("Secreted", secreted) or re.search("Secretin", secreted):
							myflag = 1
							interest_flag = 1
							mytype = unitype + "_extracellular"
							if not myid in annotation:
								annotation[myid] = []
							annotation[myid].append(mytype + "\t" + secreted + "\t" + ann + "\n" + myinfo)
				if myflag == 0:
					mytype = "NA"
					mysecreted = secreted.lower()
					if re.search("cytoplasm", mysecreted):
						mytype = unitype + "_cytoplasmic"
					if re.search("membrane", mysecreted):
						mytype = unitype + "_membrane"
					if re.search("periplasm", mysecreted):
						mytype = unitype + "_periplasmic"
					if re.search("nucleus", mysecreted):
						mytype = unitype + "_nucleus"
					if re.search("fimbrium", mysecreted):
						mytype = unitype + "_fimbrium"
					if re.search("virion", mysecreted):
						mytype = unitype + "_virion"
					if re.search("spore core", mysecreted):
						mytype = unitype + "_sporeCore"
					if re.search("mitochondrion", mysecreted):
						mytype = unitype + "_mitochondrion"
					if re.search("cell surface", mysecreted):
						mytype = unitype + "_outerMembrane"
					if re.search("cell envelop", mysecreted):
						mytype = unitype + "_cellEnvelop"
					if re.search("cell membrane", mysecreted):
						mytype = unitype + "_cellMembrane"
					if re.search("cell inner membrane", mysecreted):
						mytype = unitype + "_cellInnerMembrane"
					if re.search("cell outer membrane", mysecreted):
						mytype = unitype + "_outerMembrane"
					if re.search("cell wall", mysecreted):
						mytype = unitype + "_cellWall"
					if re.search("bacterial flagellum", mysecreted):
						mytype = unitype + "_bacterialFlagellum"
					if mytype == "NA":
						mytype = unitype + "_other"
					if mytype != "NA":
						myflag = 1
						if not myid in annotation:
							annotation[myid] = []
						annotation[myid].append(mytype + "\t" + secreted + "\t" + ann + "\n" + myinfo)

			# functional annotation
			go_bp = info[titles["GO(BP)"]]
			go_mf = info[titles["GO(MF)"]]
			go_cc = info[titles["GO(CC)"]]
			kegg = info[titles["KO"]]
			cog = info[titles["COG"]]
			pfam = info[titles["Pfam"]]
			if go_bp != "NA":
				interest_flag = 1
				mytype = unitype + "_GO(BP)"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + go_bp + "\t" + ann + "\n" + myinfo)
			if go_mf != "NA":
				interest_flag = 1
				mytype = unitype + "_GO(MF)"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + go_mf + "\t" + ann + "\n" + myinfo)
			if go_cc != "NA":
				interest_flag = 1
				mytype = unitype + "_GO(CC)"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + go_cc + "\t" + ann + "\n" + myinfo)
			if kegg != "NA":
				interest_flag = 1
				mytype = unitype + "_KEGG-KOs"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + kegg + "\t" + ann + "\n" + myinfo)
			if cog != "NA":
				interest_flag = 1
				mytype = unitype + "_COG"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + cog + "\t" + ann + "\n" + myinfo)
			if pfam != "NA":
				interest_flag = 1
				mytype = unitype + "_PfamDomain"
				if not myid in annotation:
					annotation[myid] = []
				annotation[myid].append(mytype + "\t" + pfam + "\t" + ann + "\n" + myinfo)
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
				tmp1 = "\t".join(tmp_info[titles["Length"]:len(tmp_info)])
				if member == pepid:	# use the annotation of representative
					flag = 1
					if not clust_id in outs:
						outs[clust_id] = {}
					if not clust_id in outs_info:
						outs_info[clust_id] = {}
					outs[clust_id][mytype + "\t" + mytype1 + "\t" + ann + "\t" + tmp_info[titles["Tax"]] + "\t" + tmp_info[titles["TaxID"]] + "\t" + tmp_info[titles["Rep_Tax"]] + "\t" + tmp_info[titles["Rep_TaxID"]] + "\t" + tmp_info[titles["Organism"]] + "\t" + tmp_info[titles["UniProtKB"]] + "\t" + tmp_info[titles["Entry"]] + "\t" + tmp_info[titles["Gene_names"]] + "\t" + tmp_info[titles["UniRefID"]]] = ""
					outs_info[clust_id][mytype + "\t" + myinfo] = ""
					if not clust_id in anns:
						anns[clust_id] = {}
					anns[clust_id][mytype] = ""
				# for cluster
				if not member in outs_ORF:
					outs_ORF[member] = {}
				#outs_ORF[member][myid + "\t" + tmp_info[1] + "\t" + tmp_info[2] + "\t" + tmp_info[5] + "\t" + tmp_info[4] + "\t" + tmp_info[0] + "\t" + tmp1] = ""
				outs_ORF[member][mytype + "\t" + mytype1 + "\t" + ann + "\t" + tmp_info[titles["Tax"]] + "\t" + tmp_info[titles["TaxID"]] + "\t" + tmp_info[titles["Rep_Tax"]] + "\t" + tmp_info[titles["Rep_TaxID"]] + "\t" + tmp_info[titles["Organism"]] + "\t" + tmp_info[titles["UniProtKB"]] + "\t" + tmp_info[titles["Entry"]] + "\t" + tmp_info[titles["Gene_names"]] + "\t" + tmp_info[titles["UniRefID"]]] = ""
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
				ann_info[myuniref_id][mytype + "\t" + mytype1 + "\t" + ann + "\t" + tmp_info[titles["Tax"]] + "\t" + tmp_info[titles["TaxID"]] + "\t" + tmp_info[titles["Rep_Tax"]] + "\t" + tmp_info[titles["Rep_TaxID"]] + "\t" + tmp_info[titles["Organism"]] + "\t" + tmp_info[titles["UniProtKB"]] + "\t" + tmp_info[titles["Entry"]] + "\t" + tmp_info[titles["Gene_names"]] + "\t" + tmp_info[titles["UniRefID"]]] = ""
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
				tmp2 = tmp1[titles["Length"]:len(tmp1)]
				item_num2 = len(tmp2)
				mynum = 1
				while mynum <= item_num:
					mystr = mystr + "\tNA"
					mynum = mynum + 1
				mystr2 = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
				#mynum = 1
				#while mynum <= item_num2:
				#	mystr2 = mystr2 + "\tNA"
				#	mynum = mynum + 1
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

	# get number
	number = {}
	interest = ["secreted", "cellWall", "outerMembrane", "extracellular", "signaling", "transmembrane", "GO", "KEGG-KOs", "COG"]
	other = ["cellInnerMembrane", "cellSurface", "cellEnvelop", "cellMembrane", "cytoplasm", "membrane", "cytoplasmicMembrane", "periplasm", "nucleus", "fimbrium", "virion", "sporeCore", "mitochondrion", "bacterialFlagellum"]
	unknown = ["unknown", "hypothetical"]
	unchar = ["uncharacterized"]

	total_num = 0
	for myclust in anns.keys():
		total_num = total_num + 1
		interest_flag = 0
		other_flag = 0
		unchar_flag = 0
		for mytype in anns[myclust]:
			mym, category = mytype.split("_")
			category = re.sub("GO\(BP\)", "GO", category)
			category = re.sub("GO\(MF\)", "GO", category)
			category = re.sub("GO\(CC\)", "GO", category)
			if category in interest:
				interest_flag = 1
				if not mytype in number:
					number[mytype] = {}
				number[mytype][myclust] = ""
			if category in other:
				other_flag = 1
			if category in unchar:
				unchar_flag = 1
		# foreach type of annotation
		if interest_flag == 0 and other_flag == 1:  # other type annotation
			mytype = unitype + "_others"
			if not mytype in number:
				number[mytype] = {}
			number[mytype][myclust] = ""
		if interest_flag == 0 and other_flag == 0 and unchar_flag == 1: 
			mytype = unitype + "_uncharacterized"
			if not mytype in number:
				number[mytype] = {}
			number[mytype][myclust] = ""
		if interest_flag == 0 and other_flag == 0 and unchar_flag == 0: # unknown type annotation 
			mytype = unitype + "_unknown"
			if not mytype in number:
				number[mytype] = {}
			number[mytype][myclust] = ""
    # foreach cluster

	
	#### output info ####
	outfile1 = re.sub(".detail.tsv", ".ORF.detail.tsv", outfile_detail)
	open_out = open(outfile_detail, "w")
	open_out1 = open(outfile1, "w")
	tmp1 = ann_title.split("\t")
	tmp2 = "\t".join(tmp1[titles["Length"]:len(tmp1)])
	open_out.write(utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tdescription\tTax\tTaxID\tRep_Tax\tRep_TaxID\torganism\tUniProtKB\tEntry\tGene_names\tunirefID" + "\n")
	open_out1.write(utilities.PROTEIN_ID + "\ttype\tdetail\tdescription\tTax\tTaxID\tRep_Tax\tRep_TaxID\torganism\tUniProtKB\tEntry\tGene_names\tunirefID" + "\n")
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
	
	outfile2 = re.sub(".detail.tsv", ".plot.tsv", outfile_detail)
	open_out = open(outfile2, "w")
	open_out.write("type\tnumber\tpercentage\n")
	for mytype in sorted(number.keys()):
		mynum = len(number[mytype].keys())
		myper = float(mynum)/float(total_num)
		open_out.write(mytype + "\t" + str(mynum) + "\t" + str(myper) + "\n")
	open_out.close()
	
	outfile3 = re.sub(".detail.tsv", ".mapping.tsv", outfile_detail)
	open_out = open(outfile3, "w")
	open_out.write("type\tcategory\tnumber\n")
	for mycat in sorted(maps_num.keys()):
		mytype = mycat
		mynum = len(maps_num[mycat].keys())
		open_out.write(mytype + "\t" + mycat + "\t" + str(mynum) + "\n")
	open_out.close()

	outfile4 = re.sub(".detail.tsv", ".all.spectrum.tsv", outfile_detail)
	open_file = open(outfile4, "w")
	title = "percentage\tcluster_size\tnumber"
	open_file.write(title + "\n")
	for myper in sorted(pers.keys(), key=float):
		for mytotal in sorted(pers[myper].keys(), key=int):
			mynum = 0
			mynum = len(pers[myper][mytotal].keys())
			if mynum > 0:
				open_file.write(str(myper) + "\t" + str(mytotal) + "\t" + str(mynum) + "\n")
    # foreach percentage
	open_file.close()

"""	
	outfile5 = re.sub(".detail.tsv", ".all.info.tsv", outfile_detail)
	open_file = open(outfile5, "w")
	title = "familyID\tmember\tunirefID"
	open_file.write(title + "\n")
	for myclust in sorted(outs2_info.keys()):
		for member in sorted(outs2_info[myclust].keys()):
			open_file.write(str(myclust) + "\t" + str(member) + "\t" + outs2_info[myclust][member] + "\n")
    # foreach percentage
	open_file.close()
	
	outfile6 = re.sub(".detail.tsv", ".consistency.tsv", outfile_detail)
	open_file = open(outfile6, "w")
	title = "familyID\tlabel\trep_label\tagreement"
	open_file.write(title + "\n")
	for myclust in sorted(consistency.keys()):
		open_file.write(str(myclust) + "\t" + consistency[myclust] + "\n")
    # foreach percentage
	open_file.close()
"""
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

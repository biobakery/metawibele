#!/usr/bin/env python

"""
MetaWIBELE: mspminer_protein_family module
Summary MSP annotation for each protein family

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
Summary MSP annotation for each protein family
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-m', "--msp",
	                    help='input MSP annotation file',
	                    required=True)
	parser.add_argument('-f', "--method",
	                    help='specify how to assign annotations to families',
	                    required=True,
	                    choices=["centroid", "consistency"],
	                    default="centroid")
	parser.add_argument('-o', "--output",
	                    help='output annotation detailed file for MSP annotation',
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
# assign annotation to peptide families info
#==============================================================
def collect_annotation(mspfile):	# summary_MSPminer_peptide.tsv 
	annotation = {}
	titles = {}
	open_file = open(mspfile, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	ann_title = "\t".join(info[1:len(info)])
	myindex = 0
	while myindex < len(info):
		item = info[myindex]
		titles[item] = myindex
		myindex = myindex + 1
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		if myid in annotation:
			continue
		mytype = "MSPminer_MSP"
		myinfo = "\t".join(info[1:len(info)])
		ann = info[titles["msp_name"]]
		mynote = "good"
		if ann != "msp_unknown":
			if info[-1] == "NA":
				mynote = "unclassified_MSP"
		if not myid in annotation:
			annotation[myid] = []
		annotation[myid].append(mytype + "\t" + ann + "\n" + myinfo + "\n" + mynote)
	# foreach line
	open_file.close()
	return annotation, ann_title
# collect_annotation

# assign_annotation
def assign_annotation (cutoff, pep_cluster, annotation, ann_title, assign_flag, outfile_detail):
	outs = {}
	outs_ORF = {}
	outs_info = {}
	number = {}
	titles = {}
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
				print("No MSP annotation info for the member!\t" + gene_id)
				continue
			for myall in annotation[gene_id]:
				myid, myinfo, mynote = myall.split("\n")
				mytype, ann = myid.split("\t")
				tmp_info = myinfo.split("\t")
				tmp1 = "\t".join(tmp_info[titles["gene_class"]:len(tmp_info)]) + "\t" + mynote
				if member == pepid:	# use the annotation of representative
					flag = 1
					if not clust_id in outs:
						outs[clust_id] = {}
					if not clust_id in outs_info:
						outs_info[clust_id] = {}
					outs[clust_id][myid + "\t" + tmp1] = ""
					outs_info[clust_id][myid + "\t" + myinfo] = ""
				# for cluster
				if not member in outs_ORF:
					outs_ORF[member] = {}
				outs_ORF[member][myid + "\t" + tmp1] = ""
				if tmp_info[0] != "msp_unknown":
					mytotal = mytotal + 1
				mymsp_id = tmp_info[0]
				if not clust_id in outs2_info:
					outs2_info[clust_id] = {}
				if not member in outs2_info[clust_id]:
					outs2_info[clust_id][member] = mymsp_id
				if not mymsp_id in detail_info:
					detail_info[mymsp_id] = {}
				detail_info[mymsp_id][member] = ""
				if not mymsp_id in ann_info:
					ann_info[mymsp_id] = {}
				ann_info[mymsp_id][myid + "\t" + tmp1] = ""
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
				type_tmp = "msp_unknown"
				mystr = "MSPminer_MSP"
				tmp1 = ann_title.split("\t")
				item_num = len(tmp1)
				tmp2 = tmp1[titles["gene_class"]:len(tmp1)]
				item_num2 = len(tmp2)
				mynum = 1
				while mynum <= item_num:
					mystr = mystr + "\tNA"
					mynum = mynum + 1
				mystr2 = ""
				mynum = 1
				while mynum <= item_num2:
					if mystr2 == "":
						mystr2 = "NA"
					else:
						mystr2 = mystr2 + "\tNA"
					mynum = mynum + 1
				outs[clust_id]["MSPminer_MSP" + "\t" + type_tmp + "\t" + mystr2 + "\tgood"] = ""
				outs_info[clust_id][mystr] = ""
				mymsp_id = "msp_unknown"
				for member in pep_cluster[pepid].keys():
					if not member in outs_ORF:
						outs_ORF[member] = {}
					outs_ORF[member]["MSPminer_MSP" + "\t" + type_tmp + "\t" + mystr2 + "\tgood"] = ""
					if not mymsp_id in detail_info:
						detail_info[mymsp_id] = {}
					detail_info[mymsp_id][member] = ""
				if not mymsp_id in ann_info:
					ann_info[mymsp_id] = {}
				ann_info[mymsp_id]["MSPminer_MSP" + "\t" + type_tmp + "\t" + mystr2] = ""
		# if no annotation
		cluster[clust_id] = ""

		# consistency
		detail = {}
		mytotal = len(pep_cluster[pepid].keys())
		for mytype in detail_info.keys():
			mynum = len(detail_info[mytype].keys())
		#	if mytype == "NA":
		#		continue
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
						myrep = "msp_unknown"
						if clust_id in outs:
							for myann in outs[clust_id].keys():
								tmp = myann.split("\t")
								myrep = tmp[1]	# the MSP id
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
	tmp2 = "\t".join(tmp1[titles["gene_class"]:len(tmp1)])
	open_out.write(utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\t" + tmp2 + "\tnote" + "\n")
	open_out1.write(utilities.PROTEIN_ID + "\ttype\tdetail\t" + tmp2 + "\tnote" + "\n")
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
	
	outfile4 = re.sub(".detail.tsv", ".all.spectrum.tsv", outfile_detail)
	open_file = open(outfile4, "w")
	title = "percentage\tcluster_size\tnumber"
	open_file.write(title + "\n")
	for myper in sorted(pers.keys(), key=float):
		for mytotal in sorted(pers[myper].keys(), key=int):
			mynum = 0
			for myclust in pers[myper][mytotal].keys():
				open_file.write(str(myper) + "\t" + str(mytotal) + "\t" + str(mynum) + "\n")
    # foreach percentage
	open_file.close()
	
	"""
	outfile5 = re.sub(".detail.tsv", ".all.info.tsv", outfile_detail)
	open_file = open(outfile5, "w")
	title = "familyID\tmember\tmsp_name"
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

	sys.stderr.write("### Start mspminer_protein_family.py -m " + values.msp + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	pep_cluster = collect_peptide_cluster_info (config.protein_family)
	sys.stderr.write("Get cluster info ......done\n")
	
	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	annotation, ann_title = collect_annotation(values.msp)
	sys.stderr.write("Get annotation info ......done")

	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign annotation to peptide families ......starting\n")
	assign_annotation (config.tshld_consistency, pep_cluster, annotation, ann_title, values.method, values.output)
	sys.stderr.write("\nAssign annotation to peptide families ......done\n")

	sys.stderr.write("### Finish mspminer_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

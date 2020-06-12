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
def collect_protein_cluster_info (clust_file):	# discovery_cohort.proteins.clust
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
# function collect_protein_cluster_info

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
# assign annotation to protein families info
#==============================================================
def collect_annotation(mspfile):	# summary_MSPminer_protein.tsv 
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
	cluster_component = {}
	gene_component = {}
	rep_gene = {}
	all_gene = {}
	rep_module = {}
	all_module = {}
	tmp = ann_title.split("\t")
	for item in tmp:
		titles[item] = tmp.index(item)
	for pepid in sorted(pep_cluster.keys()): # foreach protein family
		flag = 0
		detail_info = {}
		mytotal = 0
		clust_id = ""
		for member in pep_cluster[pepid].keys():
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
				mygene = tmp_info[titles["gene_class"]]
				mymodule = tmp_info[titles["module_name"]]
				if not clust_id in gene_component:
					gene_component[clust_id] = {}
				if not clust_id in rep_gene:
					rep_gene[clust_id] = {}
				if not clust_id in all_gene:
					all_gene[clust_id] = {}
				if mymodule != "NA":
					if not clust_id in all_module:
						all_module[clust_id] = {}
					all_module[clust_id][mymodule] = ""
				if re.search("core", mygene):
					if not "core" in gene_component[clust_id]:
						gene_component[clust_id]["core"] = 1
					else:
						gene_component[clust_id]["core"] = gene_component[clust_id]["core"] + 1
					if not "core" in all_gene[clust_id]:
						all_gene[clust_id]["core"] = {}
					all_gene[clust_id]["core"][member] = ""
				if re.search("accessory", mygene):
					if not "accessory" in gene_component[clust_id]:
						gene_component[clust_id]["accessory"] = 1
					else:
						gene_component[clust_id]["accessory"] = gene_component[clust_id]["accessory"] + 1
					if not "accessory" in all_gene[clust_id]:
						all_gene[clust_id]["accessory"] = {}
					all_gene[clust_id]["accessory"][member] = ""
				if member == pepid:	# use the annotation of representative
					flag = 1
					if not clust_id in outs:
						outs[clust_id] = {}
					if not clust_id in outs_info:
						outs_info[clust_id] = {}
					outs[clust_id][myid + "\t" + tmp1] = ""
					outs_info[clust_id][myid + "\t" + myinfo] = ""
					if re.search("core", mygene):
						if not "core" in rep_gene[clust_id]:
							rep_gene[clust_id]["core"] = {}
						rep_gene[clust_id]["core"][member] = ""
					if re.search("accessory", mygene):
						if not "accessory" in rep_gene[clust_id]:
							rep_gene[clust_id]["accessory"] = {}
						rep_gene[clust_id]["accessory"][member] = ""	
					if mymodule != "NA": 
						if not clust_id in rep_module:
							rep_module[clust_id] = {}
						rep_module[clust_id][mymodule] = ""
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
		tmp = {}
		mytotal = len(pep_cluster[pepid].keys())
		for mytype in detail_info.keys():
			tmp[mytype] = ""
			mynum = len(detail_info[mytype].keys())
		#	if mytype == "NA":
		#		continue
			if not mynum in detail:
				detail[mynum] = []
			detail[mynum].append(mytype)
		# foreach type
		cluster_component[clust_id] = ";".join(sorted(tmp.keys())) + "\t" + str(len(tmp.keys())) + "\t" + str(mytotal) 
		if clust_id in gene_component:
			if "core" in gene_component[clust_id]:
				gene_component[clust_id]["core"] = float(gene_component[clust_id]["core"]) / float(mytotal)
			if "accessory" in gene_component[clust_id]:
				gene_component[clust_id]["accessory"] = float(gene_component[clust_id]["accessory"]) / float(mytotal)

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
	# foreach protein cluster
	
	#### output info ####
	outfile1 = re.sub(".detail.tsv", ".ORF.detail.tsv", outfile_detail)
	outfile2 = re.sub(".detail.tsv", ".gene.detail.tsv", outfile_detail)
	outfile3 = re.sub(".detail.tsv", ".ORF.gene.detail.tsv", outfile_detail)
	open_out = open(outfile_detail, "w")
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile2, "w")
	open_out3 = open(outfile3, "w")
	tmp1 = ann_title.split("\t")
	tmp2 = "\t".join(tmp1[titles["gene_class"]:len(tmp1)])
	tmp2 = re.sub("taxa_id", "msp_component\tmulti_msp\ttotal_member\tcore_gene\taccessory_gene\tmodule_compoment\tmulti_module\ttaxa_id", tmp2)
	open_out.write(utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\t" + tmp2 + "\tnote" + "\n")
	open_out1.write(utilities.PROTEIN_ID + "\ttype\tdetail\t" + tmp2 + "\tnote" + "\n")
	open_out2.write(utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\n")
	open_out3.write(utilities.PROTEIN_ID + "\ttype\tdetail\n")
	for myclust in sorted(cluster.keys()):
		mystr = "NA\t0\t0"
		if myclust in cluster_component:
			mystr = cluster_component[myclust]
		if myclust in gene_component:
			if "core" in gene_component[myclust]:
				mystr = mystr + "\t" + str(gene_component[myclust]["core"])
			else:
				mystr = mystr + "\t0"
			if "accessory" in gene_component[myclust]:
				mystr = mystr + "\t" + str(gene_component[myclust]["accessory"])
			else:
				mystr = mystr + "\t0"
		else:
			mystr = mystr + "\t0\t0"
		if myclust in all_module:
			mytmp = sorted(all_module[myclust].keys())
			mystr = mystr + "\t" + ";".join(mytmp) + "\t" + str(len(mytmp))
		else:
			mystr = mystr + "\t" + "NA\t0" 
		tmp = mystr.split("\t")
		multi = tmp[0].split(";")
		if len(multi) > 1:
			multi = tmp[0]
		else:
			multi = "NA"
		rep_core = "NA"
		rep_acc = "NA"
		module1 = "NA"
		all_core = "NA"
		all_acc = "NA"
		module2 = "NA"
		if myclust in rep_gene:
			if "core" in rep_gene[myclust]:
				rep_core = ";".join(sorted(rep_gene[myclust]["core"].keys()))
			if "accessory" in rep_gene[myclust]:
				rep_acc = ";".join(sorted(rep_gene[myclust]["accessory"].keys()))
		if myclust in rep_module:
			module1 = ";".join(sorted(rep_module[myclust].keys()))
		if myclust in all_gene:
			if "core" in all_gene[myclust]:
				all_core = ";".join(sorted(all_gene[myclust]["core"].keys()))
			if "accessory" in all_gene[myclust]:
				all_acc = ";".join(sorted(all_gene[myclust]["accessory"].keys()))
		if myclust in all_module:
			module2 = ";".join(sorted(all_module[myclust].keys()))
		if multi != "NA":
			open_out2.write(myclust + "\t" + "MSPminer_multi-msp" + "\t" + multi + "\n")
		if rep_core != "NA":
			open_out2.write(myclust + "\t" + "MSPminer_core" + "\t" + rep_core + "\n")
		if all_core != "NA":
			open_out2.write(myclust + "\t" + "MSPminer_family-core" + "\t" + all_core + "\n")
		if rep_acc != "NA":
			open_out2.write(myclust + "\t" + "MSPminer_accessory" + "\t" + rep_acc + "\n")
		if all_acc != "NA":
			open_out2.write(myclust + "\t" + "MSPminer_family-accessory" + "\t" + all_acc + "\n")
		if module1 != "NA":
			open_out2.write(myclust + "\t" + "MSPminer_module" + "\t" + module1 + "\n")
		if module2 != "NA":
			open_out2.write(myclust + "\t" + "MSPminer_family-module" + "\t" + module2 + "\n")
			mytmp = module2.split(";")
			if len(mytmp) > 1:
				open_out2.write(myclust + "\t" + "MSPminer_multi-module" + "\t" + module2 + "\n")

		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in outs:
				for myid in sorted(outs[myclust].keys()):
					mytmp = myid.split("\t")
					mytmp[4] = mytmp[4] + "\t" + mystr
					myid = "\t".join(mytmp)
					open_out.write(myclust + "\t" + myid + "\n")
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in outs2:
				for myid in sorted(outs2[myclust].keys()):
					mytmp = myid.split("\t")
					mytmp[4] = mytmp[4] + "\t" + mystr
					myid = "\t".join(mytmp)
					open_out.write(myclust + "\t" + myid + "\n")
	# foreach cluster
	open_out.close()
	open_out2.close()
	for mypep in sorted(outs_ORF.keys()):
		for myid in sorted(outs_ORF[mypep].keys()):
			mytmp = myid.split("\t") 
			mygene = mytmp[2]
			mymodule = mytmp[3]
			mytmp[4] = mytmp[4] + "\t" + mytmp[1] + "\t1\t1"
			if re.search("core", mygene):
				open_out3.write(mypep + "\t" + "MSPminer_core" + "\t" + mypep + "\n")
				mytmp[4] = mytmp[4] + "\t" + "1\t0"	
			if re.search("accessory", mygene):
				open_out3.write(mypep + "\t" + "MSPminer_accessory" + "\t" + mypep + "\n")
				mytmp[4] = mytmp[4] + "\t" + "0\t1"
			if mygene == "NA":
				mytmp[4] = mytmp[4] + "\t" + "0\t0"
			if mymodule != "NA":
				open_out3.write(mypep + "\t" + "MSPminer_module" + "\t" + mymodule + "\n")
			mytmp[4] = mytmp[4] + "\t" + mymodule + "\t1"
			myid = "\t".join(mytmp)
			open_out1.write(mypep + "\t" + myid + "\n")
		# foreach type
	# foreach protein family
	open_out1.close()
	open_out3.close()

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
	pep_cluster = collect_protein_cluster_info (config.protein_family)
	sys.stderr.write("Get cluster info ......done\n")
	
	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	annotation, ann_title = collect_annotation(values.msp)
	sys.stderr.write("Get annotation info ......done")

	### assign annotation to protein families ###
	sys.stderr.write("\nAssign annotation to protein families ......starting\n")
	assign_annotation (config.tshld_consistency, pep_cluster, annotation, ann_title, values.method, values.output)
	sys.stderr.write("\nAssign annotation to protein families ......done\n")

	sys.stderr.write("### Finish mspminer_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python

"""
MetaWIBELE: psortb_protein_family module
Summary the results of Psortb with Gram+, Gram- mode, Archaea mode

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
Summary the results of Psortb with Gram+, Gram- mode, Archaea mode
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
	                    default="psortb.gram_positive.out.location.tsv")
	parser.add_argument('-p', "--path",
	                    help='input the path of annotation file',
	                    required=True)
	parser.add_argument('-a', "--method",
	                    help='specify how to assign annotations to families',
	                    choices=["centroid", "consistency"],
	                    required=True,
	                    default="consistency")
	parser.add_argument('-o', "--output",
	                    help='output annotation summary file',
	                    required=True)
	values = parser.parse_args()

	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_cluster_info (clust_file):  # discovery_cohort.peptides.clust
	cluster = {}
	cluster_mem = {}
	open_file = open(clust_file, "r")
	myclust = ""
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
		cluster_mem[myid] = myclust + "\t" + myclust_id
		cluster[myclust][myid] = myclust_id
    # foreach line
	open_file.close()
	return cluster, cluster_mem
# function collect_cluster_info


#==============================================================
# collect localizing info
#==============================================================
def collect_localizing_info (cluster_mem, extension, ann_path, outfile):	# list.txt
	gram_p = {}
	gram_n = {}
	archaea = {}
	location = {}
	score = {}
	location_n = {}
	location_p = {}
	location_a = {}
	filelist = utilities.find_files(ann_path, extension, None)
	for myfile in filelist:
		#myfile = ann_path + "/" + samplelist + "/" + samplelist + ".psortb.gram_positive.out.location.tsv"
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			open_file = open(myfile, "r")
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^name", line): 
					continue
				info = line.split("\t")
				myid = info[0]
				if not myid in cluster_mem:
					continue
				sample = re.sub("_[\d]+$", "", myid)
				if not sample in gram_p:
					gram_p[sample] = {}
				if not info[0] in location:
					location[info[0]] = info[1]
					score[info[0]] = info[2]
				else:
					if float(info[2]) > float(score[info[0]]):
						location[info[0]] = info[1]
						score[info[0]] = info[2]
				location_p[info[0]] = info[1]
				gram_p[sample][info[0]] = info[1]
			# foreach line
			open_file.close()
			
		#myfile = ann_path + "/" + samplelist + "/" + samplelist + ".psortb.gram_negative.out.location.tsv"
		myfile1 = re.sub("psortb.gram_positive.out.location.tsv", "psortb.gram_negative.out.location.tsv", myfile)
		if not os.path.isfile(myfile1):
			print ("File not exist!\t" + myfile1)
		else:
			open_file = open(myfile1, "r")
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^name", line): 
					continue
				info = line.split("\t")
				myid = info[0]
				if not myid in cluster_mem:
					continue
				sample = re.sub("_[\d]+$", "", myid)
				if not sample in gram_n:
					gram_n[sample] = {}
				if not info[0] in location:
					location[info[0]] = info[1]
					score[info[0]] = info[2]
				else:
					if float(info[2]) > float(score[info[0]]):
						location[info[0]] = info[1]
						score[info[0]] = info[2]
				location_n[info[0]] = info[1]
				gram_n[sample][info[0]] = info[1]
			# foreach line
			open_file.close()

		#myfile = ann_path + "/" + samplelist + "/" + samplelist + ".psortb.archaea.out.location.tsv"
		myfile1 = re.sub("psortb.gram_positive.out.location.tsv", "psortb.archaea.out.location.tsv", myfile)
		if not os.path.isfile(myfile1):
			print ("File not exist!\t" + myfile1)
		else:
			open_file = open(myfile1, "r")
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^name", line): 
					continue
				info = line.split("\t")
				myid = info[0]
				if not myid in cluster_mem:
					continue
				sample = re.sub("_[\d]+$", "", myid)
				if not sample in archaea:
					archaea[sample] = {}
				if not info[0] in location:
					location[info[0]] = info[1]
					score[info[0]] = info[2]
				else:
					if float(info[2]) > float(score[info[0]]):
						location[info[0]] = info[1]
						score[info[0]] = info[2]
				location_a[info[0]] = info[1]
				archaea[sample][info[0]] = info[1]
			# foreach line
			open_file.close()
	# foreach samplelist
	
	# output details
	outfile1 = re.sub("_proteinfamilies.", "_proteinfamilies.ORF.", outfile)
	open_out = open(outfile1, "w")
	open_out.write(utilities.PROTEIN_ID + "\ttype\tdetail\tscore\n")
	for myid in sorted(location.keys()):
		mytype = location[myid]
		if mytype == "Unknown":
			mytype = "PSORTb_unknown"
		if mytype == "Cytoplasmic":
			mytype = "PSORTb_cytoplasmic"
		if mytype == "CytoplasmicMembrane":
			mytype = "PSORTb_cytoplasmicMembrane"
		if mytype == "Extracellular":
			mytype = "PSORTb_extracellular"
		if mytype == "Cellwall":
			mytype = "PSORTb_cellWall"
		if mytype == "OuterMembrane":
			mytype = "PSORTb_outerMembrane"
		if mytype == "Periplasmic":
			mytype = "PSORTb_periplasmic"
		myscore = "NA"
		if myid in score:
			myscore = score[myid]
		mydetail = re.sub("PSORTb_", "", mytype)
		open_out.write(myid + "\t" + mytype + "\t" + mydetail + "\t" + str(myscore) + "\n")
	# foreach seqID
	open_out.close()

	return gram_p, gram_n, archaea, location, score, location_p, location_n, location_a
# function collect_localizing_info


#==============================================================
# output info
#==============================================================
def output_info (cutoff, cluster, gram_p, gram_n, archaea, location, score,location_p, location_n, location_a, assign_flag, outfile):

	"""
	open_file = open(outfile, "w")
	title = "sample\tdiagnosis\tmode\tseqID\ttype\tnumber"
	open_file.write(title + "\n")
	for sample in sorted(gram_p.keys()):
		mymode = "Gram+"
		mydia = "NA"
		if sample in samples:
			mydia = samples[sample]
		for myid in sorted(gram_p[sample].keys()):
			mystr = sample + "\t" + mydia + "\t" + mymode + "\t" + myid + "\t" + gram_p[sample][myid]
			open_file.write(mystr + "\t1\n")
	# foreach sample
	
	for sample in sorted(gram_n.keys()):
		mymode = "Gram-"
		mydia = "NA"
		if sample in samples:
			mydia = samples[sample]
		for myid in sorted(gram_n[sample].keys()):
			mystr = sample + "\t" + mydia + "\t" + mymode + "\t" + myid + "\t" + gram_n[sample][myid]
			open_file.write(mystr + "\t1\n")
	# foreach sample
	
	for sample in sorted(archaea.keys()):
		mymode = "Archaea"
		mydia = "NA"
		if sample in samples:
			mydia = samples[sample]
		for myid in sorted(archaea[sample].keys()):
			mystr = sample + "\t" + mydia + "\t" + mymode + "\t" + myid + "\t" + archaea[sample][myid]
			open_file.write(mystr + "\t1\n")
	# foreach sample
	open_file.close()

	outfile2_1 = re.sub(".tsv", ".member_detail.tsv", outfile)
	open_file2_1 = open(outfile2_1, "w")
	open_file2_1.write("familyID\tmember\ttype\n")
	"""

	open_file2 = open(outfile, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tconsistency"
	open_file2.write(title + "\n")
	pers = {}
	ann_info = {}
	for myclust in cluster.keys():
		mytotal = len(cluster[myclust].keys())
		detail_info = {}
		detail = {}
		detail_rep = {}
		detail_info_rep = {}
		clust_id = ""
		for myid in cluster[myclust].keys():
			clust_id = cluster[myclust][myid]	
			if not clust_id in ann_info:
				ann_info[clust_id] = {}
			if myid in location:
				mytype = location[myid]
				if mytype == "Unknown":
					mytype = "PSORTb_unknown"
				if mytype == "Cytoplasmic":
					mytype = "PSORTb_cytoplasmic"
				if mytype == "CytoplasmicMembrane":
					mytype = "PSORTb_cytoplasmicMembrane"
				if mytype == "Extracellular":
					mytype = "PSORTb_extracellular"
				if mytype == "Cellwall":
					mytype = "PSORTb_cellWall"
				if mytype == "OuterMembrane":
					mytype = "PSORTb_outerMembrane"
				if mytype == "Periplasmic":
					mytype = "PSORTb_periplasmic"
				if not mytype in detail_info:
					detail_info[mytype] = {}
				detail_info[mytype][myid] = ""
				if myclust == myid:
					if not mytype in detail_info_rep:
						detail_info_rep[mytype] = {}
					detail_info_rep[mytype][myid] = ""
				if not myid in ann_info[clust_id]:
					ann_info[clust_id][myid] = {}
				ann_info[clust_id][myid][mytype] = ""
			# if labeling myid
			else:
				if not myid in ann_info[clust_id]:
					ann_info[clust_id][myid] = {}
				ann_info[clust_id][myid]["PSORTb_no"] = ""
		# froeach member

		for mytype in detail_info.keys():
			mynum = len(detail_info[mytype].keys())
			if not mynum in detail:
				detail[mynum] = []
			detail[mynum].append(mytype)
		for mytype in detail_info_rep.keys():
			mynum = len(detail_info_rep[mytype].keys())
			if not mynum in detail_rep:
				detail_rep[mynum] = []
			detail_rep[mynum].append(mytype)
		# foreach type
		if assign_flag == "centroid":  # assign annotation based on representative info
			for mynum in sorted(detail_rep.keys(), key=int, reverse=True):
				item_num = len(detail_rep[mynum])
				for item in detail_rep[mynum]:
					if item_num >=2 and re.search("_unknown", item):
						continue
					if item == "PSORTb_cytoplasmic" and "PSORTb_cytoplasmicMembrane" in detail_rep[mynum]:
						continue
					myper = round(float(mynum)/float(mytotal), 2)
					if mynum != 0:
						if not myper in pers:
							pers[myper] = {}
						if not mytotal in pers[myper]:
							pers[myper][mytotal] = {}
						pers[myper][mytotal][clust_id] = ""
					mydetail = re.sub("PSORTb_", "", item)
					#if myper >= float(cutoff):
					open_file2.write(clust_id + "\t" + item + "\t" + mydetail + "\t" + str(myper) + "\n")
				# foreach type
				break
			# foreach voting weight
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			for mynum in sorted(detail.keys(), key=int, reverse=True):
				item_num = len(detail[mynum])
				for item in detail[mynum]:
					if item_num >=2 and re.search("_unknown", item):
						continue
					if item == "PSORTb_cytoplasmic" and "PSORTb_cytoplasmicMembrane" in detail[mynum]:
						continue
					myper = round(float(mynum)/float(mytotal), 2)
					if mynum != 0:
						if not myper in pers:
							pers[myper] = {}
						if not mytotal in pers[myper]:
							pers[myper][mytotal] = {}
						pers[myper][mytotal][clust_id] = ""
					if myper >= float(cutoff):
						mydetail = re.sub("PSORTb_", "", item)
						open_file2.write(clust_id + "\t" + item + "\t" + mydetail + "\t" + str(myper) + "\n")
				# foreach type
				break
			# foreach voting weight
		# if consistency
	# foreach clusters
	open_file2.close()

	outfile2 = re.sub(".detail.tsv", ".all.spectrum.tsv", outfile)
	open_file2 = open(outfile2, "w")
	title = "percentage\tmember_num\tnumber"
	open_file2.write(title + "\n")
	for myper in sorted(pers.keys(), key=float):
		for mytotal in sorted(pers[myper].keys(), key=int):
			mynum = len(pers[myper][mytotal].keys())
			open_file2.write(str(myper) + "\t" + str(mytotal) + "\t" + str(mynum) + "\n")
	# foreach percentage
	open_file2.close()

	"""
	outfile2 = re.sub(".tsv", ".detail.plot.tsv", outfile)
	open_file2 = open(outfile2, "w")
	title = "Gram_p\tGram_n\tArchaea"
	open_file2.write(title + "\n")
	for myid in sorted(location_p.keys()):
		mystr = myid
		if myid in location_n:
			mystr = mystr + "\t" + myid
		else:
			mystr = mystr + "\tNA"
		if myid in location_a:
			mystr = mystr + "\t" + myid
		else:
			mystr = mystr + "\tNA"
		open_file2.write(mystr + "\n")
	for myid in sorted(location_n.keys()):
		if myid in location_p:
			continue
		else:
			if myid in location_a:
				mystr = "NA" + "\t" + myid + "\t" + myid
			else:
				mystr = "NA" + "\t" + myid + "\tNA"
			open_file2.write(mystr + "\n")
	for myid in sorted(location_a.keys()):
		if myid in location_p:
			continue
		if myid in location_n:
			continue
		mystr = "NA" + "\tNA\t" + myid
		open_file2.write(mystr + "\n")
	open_file2.close()
	"""

# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start psortb_protein_family.py ####\n")
	
	### collect sample info ###
	sys.stderr.write("Get info ......starting\n")
	cluster, cluster_mem = collect_cluster_info (config.protein_family)
	sys.stderr.write("Get info ......done\n")

	### collect annotation info and do stat ###
	sys.stderr.write("Get localizing info ......starting\n")
	gram_p, gram_n, archaea, location, scores, location_p, location_n, location_a = collect_localizing_info (cluster_mem, values.extension, values.path, values.output)
	sys.stderr.write("Get localizing info ......done\n")
	
	### Output sample info
	sys.stderr.write("\nOutput localizing summary info ......starting\n")
	output_info (config.tshld_consistency, cluster, gram_p, gram_n, archaea, location, scores, location_p, location_n, location_a, values.method, values.output)
	sys.stderr.write("Output localizing summary info ......done\n")

	sys.stderr.write("### Finish psortb_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

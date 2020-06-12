#!/usr/bin/env python

"""
MetaWIBELE: interproscan_tmhmm_protein_family module
Summary the results of TMHMM from InterProScan

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
Summary the results of TMHMM from InterProScan
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
	                    default="tmhmm.transmembrane.tsv")
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
# collect transmembrane info
#==============================================================
def collect_transmembrane_info (cluster_mem, extension, ann_path, outfile): # list.txt
	transmem = {}
	details = {}
	filelist = utilities.find_files(ann_path, extension, None)
	for myfile in filelist:
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			open_file = open(myfile, "r")
			titles = {}
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				info = line.split("\t")
				if re.search("^" + utilities.PROTEIN_ID, line):
					for item in info:
						titles[item] = info.index(item)
					continue
				myid = info[titles[utilities.PROTEIN_ID]]
				if not myid in cluster_mem:
					continue
				details[info[titles[utilities.PROTEIN_ID]]] = info[titles["Prediction"]]
				sample = re.sub("_[\d]+$", "", myid)
				if not sample in transmem:
					transmem[sample] = {}
				transmem[sample][info[titles[utilities.PROTEIN_ID]]] = info[titles["Prediction"]]
			# foreach line
			open_file.close()
	# foreach samplelist

	# output details
	outfile1 = re.sub("_proteinfamilies.", "_proteinfamilies.ORF.", outfile)
	open_out = open(outfile1, "w")
	open_out.write(utilities.PROTEIN_ID + "\ttype\tdetail\tdescription\n")
	for myid in sorted(details.keys()):
		myinfo = details[myid]
		open_out.write(myid + "\tTMHMM_transmembrane\tTMHMM_transmembrane\t" + myinfo + "\n")
	open_out.close()
	
	return transmem
# function collect_transmembrane_info


#==============================================================
# output info
#==============================================================
def output_info (cutoff, cluster, transmem, assign_flag, outfile):
	annotation = {}
	for sample in sorted(transmem.keys()):
		for myid in sorted(transmem[sample].keys()):
			annotation[myid] = "TMHMM_transmembrane" + "\t" + str(transmem[sample][myid])
			#mystr = sample + "\t" + mydia + "\t" + myid + "\t" + str(transmem[sample][myid])
			#open_file.write(mystr + "\t1\n")
	# foreach sample

	open_file2 = open(outfile, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tconsistency"
	details = {}
	details_rep = {}
	cluster_id = {}
	pers = {}
	pers_cluster = {}
	ann_info = {}
	open_file2.write(title + "\n")
	for myclust in cluster.keys():
		mytotal = len(cluster[myclust].keys())
		mynum = 0
		clust_id = ""
		for myid in cluster[myclust].keys():
			clust_id = cluster[myclust][myid]
			if not clust_id in ann_info:
				ann_info[clust_id] = {}
			if myid in annotation:
				mynum = mynum + 1
				ann_info[clust_id][myid] = "TMHMM_transmembrane"
			else:
				ann_info[clust_id][myid] = "No"
		# foreach member
		myper = round(float(mynum)/float(mytotal), 2)
		if mynum != 0:
			cluster_id[clust_id] = ""
			if not myper in pers:
				pers[myper] = {}
			if not mytotal in pers[myper]:
				pers[myper][mytotal] = {}
			pers[myper][mytotal][clust_id] = ""
			pers_cluster[clust_id] = str(mytotal) + "\t" + str(mynum)
			if myper >= float(cutoff):
				details[clust_id] = "TMHMM_transmembrane\tTMHMM_transmembrane" + "\t" + str(myper)
			if myclust in annotation:  # the centroid is annotated
				details_rep[clust_id] = "TMHMM_transmembrane\tTMHMM_transmembrane" + "\t" + str(myper)
	# foreach cluster
	for myclust in sorted(cluster_id.keys()):
		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in details_rep:
				open_file2.write(myclust + "\t" + details_rep[myclust] + "\n")
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in details:
				open_file2.write(myclust + "\t" + details[myclust] + "\n")
    # foreach cluster
	#open_file2.close()
	#for myclust in sorted(ann_info.keys()):
	#	for myid in sorted(ann_info[myclust].keys()):
	#		open_file2_1.write(myclust + "\t" + myid + "\t" + ann_info[myclust][myid] + "\n")
	#open_file2_1.close()
	
	outfile2 = re.sub(".detail.tsv", ".all.spectrum.tsv", outfile)
	open_file2 = open(outfile2, "w")
	title = "percentage\tmember_num\tnumber"
	open_file2.write(title + "\n")
	for myper in sorted(pers.keys(), key=float):
		for mytotal in sorted(pers[myper].keys(), key=int):
			mynum = len(pers[myper][mytotal].keys())
			open_file2.write(str(myper) + "\t" + str(mytotal) + "\t" + str(mynum) + "\n")
	open_file2.close()

	"""
	outfile2 = re.sub(".tsv", ".all.cluster_info.tsv", outfile)
	open_file2 = open(outfile2, "w")
	title = "familyID\tmember_num\tmotif_num"
	open_file2.write(title + "\n")
	for myid in sorted(pers_cluster.keys()):
		open_file2.write(myid + "\t" + pers_cluster[myid] + "\n")
	open_file2.close()
	"""

# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start interproscan_tmhmm_protein_family.py -p " + values.path + " ####\n")
	
	### collect sample info ###
	sys.stderr.write("Get info ......starting\n")
	cluster, cluster_mem = collect_cluster_info (config.protein_family)
	sys.stderr.write("Get info ......done\n")

	### collect annotation info and do stat ###
	sys.stderr.write("Get transmembrane info ......starting\n")
	transmem = collect_transmembrane_info (cluster_mem, values.extension, values.path, values.output)
	sys.stderr.write("Get transmembrane info ......done\n")
	
	### Output sample info
	sys.stderr.write("\nOutput transmembrane summary info ......starting\n")
	output_info (config.tshld_consistency, cluster, transmem, values.method, values.output)
	sys.stderr.write("Output transmembrane summary info ......done\n")

	sys.stderr.write("### Finish interproscan_tmhmm_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

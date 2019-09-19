#!/usr/bin/env python

"""
MetaWIBELE: interproscan_phobius_protein_family module
Summary the results of Phobius from InterProScan

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
	#import config
	#import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Summary the results of Phobius from InterProScan
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
	                    default="phobius.signaling.tsv")
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
# collect phobius info
#==============================================================
def collect_phobius_info (cluster_mem, extension, ann_path, outfile):	# list.txt
	transmem = {}
	signal = {}
	detail_signal = {}
	detail_trans = {}
	filelist = utilities.find_files(ann_path, extension, None)
	for myfile in filelist:
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			open_file = open(myfile, "r")
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^" + utilities.PROTEIN_ID, line):
					continue
				info = line.split("\t")
				myid = info[0]
				if not myid in cluster_mem:
					continue
				sample = re.sub("_[\d]+$", "", myid)
				if not sample in signal:
					signal[sample] = {}
				signal[sample][info[0]] = info[2]
				detail_signal[info[0]] = info[2]
			# foreach line
			open_file.close()

		myfile = re.sub(extension, "phobius.transmembrane.tsv", myfile)
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			open_file = open(myfile, "r")
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^" + utilities.PROTEIN_ID, line):
					continue
				info = line.split("\t")
				myid = info[0]
				if not myid in cluster_mem:
					continue
				sample = re.sub("_[\d]+$", "", myid)
				if not sample in transmem:
					transmem[sample] = {}
				transmem[sample][info[0]] = info[2]
				detail_trans[info[0]] = info[2]
			# foreach line
	# foreach samplelist

	# output details
	outfile1 = re.sub("_proteinfamilies.", "_proteinfamilies.ORF.signaling.", outfile)
	open_out = open(outfile1, "w")
	open_out.write(utilities.PROTEIN_ID + "\ttype\tdetail\tdescription\n")
	for myid in sorted(detail_signal.keys()):
		myinfo = detail_signal[myid]
		open_out.write(myid + "\tPhobius_signaling\tPhobius_signaling\t" + myinfo + "\n")
	# foreach seqID
	open_out.close()
	outfile1 = re.sub("_proteinfamilies.", "_proteinfamilies.ORF.transmembrane.", outfile)
	open_out = open(outfile1, "w")
	open_out.write(utilities.PROTEIN_ID + "\ttype\tdetail\tdescription\n")
	for myid in sorted(detail_trans.keys()):
		myinfo = detail_trans[myid]
		open_out.write(myid + "\tPhobius_transmembrane\tPhobius_transmembrane\t" + myinfo + "\n")
	# foreach seqID
	open_out.close()

	return signal, transmem
# function collect_phobius_info


#==============================================================
# output info
#==============================================================
def output_info (cutoff_s, cutoff_t, cluster, signal, transmem, assign_flag, outfile):

	"""
	outfile1 = re.sub(".tsv", ".signaling.tsv", outfile)
	open_file = open(outfile1, "w")
	title = "sample\tdiagnosis\tseqID\tSP\tNumber"
	open_file.write(title + "\n")
	"""
	annotation_s = {}
	for sample in sorted(signal.keys()):
		for myid in sorted(signal[sample].keys()):
			annotation_s[myid] = "Phobius_signaling" + "\t" + str(signal[sample][myid])
			#mystr = sample + "\t" + mydia + "\t" + myid + "\t" + str(signal[sample][myid])
			#open_file.write(mystr + "\t1\n")
	# foreach sample
	
	"""
	outfile2 = re.sub(".tsv", ".transmembrane.tsv", outfile)
	open_file = open(outfile2, "w")
	title = "sample\tdiagnosis\tseqID\tTM\tNumber"
	open_file.write(title + "\n")
	"""
	annotation_t = {}
	for sample in sorted(transmem.keys()):
		for myid in sorted(transmem[sample].keys()):
			annotation_t[myid] = "Phobius_transmembrane" + "\t" + str(transmem[sample][myid])
			#mystr = sample + "\t" + mydia + "\t" + myid + "\t" + str(transmem[sample][myid])
			#open_file.write(mystr + "\t1\n")
	# foreach sample

	# get details and numbers
	detail_s = {}
	detail_s_rep = {}
	cluster_id_s = {}
	pers_s = {}
	pers_s_cluster = {}
	detail_t = {}
	detail_t_rep = {}
	cluster_id_t = {}
	pers_t = {}
	pers_t_cluster = {}
	for myclust in cluster.keys():
		mytotal = len(cluster[myclust].keys())
		mynum_s = 0
		mynum_t = 0
		clust_id = ""
		for myid in cluster[myclust].keys():
			clust_id = cluster[myclust][myid]
			if myid in annotation_s:
				mynum_s = mynum_s + 1
			if myid in annotation_t:
				mynum_t = mynum_t + 1
		# foreach member
		myper_s = round(float(mynum_s)/float(mytotal), 2)
		myper_t = round(float(mynum_t)/float(mytotal), 2)
		if mynum_s != 0:
			cluster_id_s[clust_id] = ""
			if not myper_s in pers_s:
				pers_s[myper_s] = {}
			if not mytotal in pers_s[myper_s]:
				pers_s[myper_s][mytotal] = {}
			pers_s[myper_s][mytotal][clust_id] = ""
			pers_s_cluster[clust_id] = str(mytotal) + "\t" + str(mynum_s)
			if myper_s >= float(cutoff_s):
				detail_s[clust_id] = "Phobius_signaling\tPhobius_signaling" + "\t" + str(myper_s)
			if myclust in annotation_s:  # the centroid is annotated
				detail_s_rep[clust_id] = "Phobius_signaling\tPhobius_signaling" + "\t" + str(myper_s)
		if mynum_t != 0:
			cluster_id_t[clust_id] = ""
			if not myper_t in pers_t:
				pers_t[myper_t] = {}
			if not mytotal in pers_t[myper_t]:
				pers_t[myper_t][mytotal] = {}
			pers_t[myper_t][mytotal][clust_id] = ""
			pers_t_cluster[clust_id] = str(mytotal) + "\t" + str(mynum_t)
			if myper_t >= float(cutoff_t):
				detail_t[clust_id] = "Phobius_transmembrane\tPhobius_transmembrane" + "\t" + str(myper_t)
			if myclust in annotation_t:  # the centroid is annotated
				detail_t_rep[clust_id] = "Phobius_transmembrane\tPhobius_transmembrane" + "\t" + str(myper_t)
	# foreach cluster

	outfile3 = re.sub(".detail.tsv", ".signaling.detail.tsv", outfile)
	open_file3 = open(outfile3, "w")
	outfile4 = re.sub(".detail.tsv", ".transmembrane.detail.tsv", outfile)
	open_file4 = open(outfile4, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tconsistency"
	open_file3.write(title + "\n")
	open_file4.write(title + "\n")
	for myclust in sorted(cluster_id_s.keys()):
		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in detail_s_rep:
				open_file3.write(myclust + "\t" + detail_s_rep[myclust] + "\n")
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in detail_s:
				open_file3.write(myclust + "\t" + detail_s[myclust] + "\n")
    # foreach cluster
	open_file3.close()
	for myclust in sorted(cluster_id_t.keys()):
		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in detail_t_rep:
				open_file4.write(myclust + "\t" + detail_t_rep[myclust] + "\n")
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in detail_t:
				open_file4.write(myclust + "\t" + detail_t[myclust] + "\n")
	# foreach cluster
	open_file4.close()
	
	outfile5 = re.sub(".detail.tsv", ".signaling.all.spectrum.tsv", outfile)
	open_file5 = open(outfile5, "w")
	outfile6 = re.sub(".detail.tsv", ".transmembrane.all.spectrum.tsv", outfile)
	open_file6 = open(outfile6, "w")
	title = "percentage\tmember_num\tnumber"
	open_file5.write(title + "\n")
	open_file6.write(title + "\n")
	for myper in sorted(pers_s.keys(), key=float):
		for mytotal in sorted(pers_s[myper].keys(), key=int):
			mynum = len(pers_s[myper][mytotal].keys())
			open_file5.write(str(myper) + "\t" + str(mytotal) + "\t" + str(mynum) + "\n")
	open_file5.close()
	for myper in sorted(pers_t.keys(), key=float):
		for mytotal in sorted(pers_t[myper].keys(), key=int):
			mynum = len(pers_t[myper][mytotal].keys())
			open_file6.write(str(myper) + "\t" + str(mytotal) + "\t" + str(mynum) + "\n")
	open_file6.close()
	
	"""
	outfile5 = re.sub(".tsv", ".signaling.all.cluster_info.tsv", outfile)
	open_file5 = open(outfile5, "w")
	outfile6 = re.sub(".tsv", ".transmembrane.all.cluster_info.tsv", outfile)
	open_file6 = open(outfile6, "w")
	title = "familyID\tmember_num\tmotif_num"
	open_file5.write(title + "\n")
	open_file6.write(title + "\n")
	for myid in sorted(pers_s_cluster.keys()):
		open_file5.write(myid + "\t" + pers_s_cluster[myid] + "\n")
	open_file5.close()
	for myid in sorted(pers_t_cluster.keys()):
		open_file6.write(myid + "\t" + pers_t_cluster[myid] + "\n")
	open_file6.close()
	"""

# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start interproscan_phobius_protein_family.py -p " + values.path + " ####\n")
	
	### collect sample info ###
	sys.stderr.write("Get info ......starting\n")
	cluster, cluster_mem = collect_cluster_info (config.protein_family)
	sys.stderr.write("Get info ......done\n")

	### collect annotation info and do stat ###
	sys.stderr.write("Get phobius info ......starting\n")
	signal, transmem = collect_phobius_info (cluster_mem, values.extension, values.path, values.output)
	sys.stderr.write("Get phobius info ......done\n")
	
	### Output sample info
	sys.stderr.write("\nOutput signaling and transmembrane summary info ......starting\n")
	output_info (config.tshld_consistency, config.tshld_consistency, cluster, signal, transmem, values.method, values.output)
	sys.stderr.write("Output signaling and transmembrane summary info ......done\n")

	sys.stderr.write("### Finish interproscan_phobius_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

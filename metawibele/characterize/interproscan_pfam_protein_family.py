#!/usr/bin/env python

"""
MetaWIBELE: interproscan_pfam_protein_family module
Summary the results of Pfam annotation from InterProScan for families
Configuration settings

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
Summary the results of Pfam annotation from InterProScan
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-e', "--extension",
	                    help='the file extension to search for',
	                    required=True,
						default="interpro.PfamDomain.tsv")
	parser.add_argument('-p', "--path",
	                    help='the path of folder',
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
# collect Pfam info
#==============================================================
def collect_pfam_ann (pfamfile):
	pfam_ann = {}
	for line in utils.gzip_bzip2_biom_open_readlines (pfamfile): 
		line = line.strip()
		if not len(line):
			continue
		if re.search("^Pfam", line):
			continue
		info = line.split("\t")
		pfam_ann[info[0]] = info[1]
	# foreach line
	
	return pfam_ann
# collect_pfam_ann


#==============================================================
# collect cluster info
#==============================================================
def collect_cluster_info (clust_file):  # discovery_cohort.peptides.clust
	cluster = {}
	cluster_id = {}
	cluster_mem = {}
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
			cluster_id[myclust] = myclust_id
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster_mem[myid] = myclust + "\t" + myclust_id
		cluster[myclust][myid] = myclust_id
    # foreach line
	open_file.close()
	return cluster, cluster_id, cluster_mem
# function collect_cluster_info


#==============================================================
# collect pfam info
#==============================================================
def collect_pfam_info (cluster_mem, extension, ann_path, outfile):	# list.txt
	pfams = {}
	filelist = utilities.find_files(ann_path, extension, None)
	for myfile in filelist:
		#myfile = ann_path + "/" + samplelist + "/" + samplelist + ".interpro.PfamDomain.tsv"
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
			continue
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
			pfam = info[1]
			if not myid in pfams:
				pfams[myid] = {}
			pfams[myid][pfam] = info[2]
		# foreach line
		open_file.close()
	# foreach samplelist

	# output details
	outfile1 = re.sub("_proteinfamilies.", "_proteinfamilies.ORF.", outfile)
	open_out = open(outfile1, "w")
	open_out.write(utilities.PROTEIN_ID + "\ttype\tdetail\tdescription\n")
	for myid in sorted(pfams.keys()):
		myacc = ""
		myann = ""
		for item in sorted(pfams[myid].keys()):
			if myacc == "":
				myacc = item
			else:
				myacc = myacc + ";" + item
			if myann == "":
				myann = pfams[myid][item]
			else:
				myann = myann + ";" + pfams[myid][item]
		open_out.write(myid + "\tPfam_PfamDomain\t" + myacc + "\t" + myann + "\n")
	# foreach seqID
	open_out.close()

	return pfams
# function collect_pfam_info


#==============================================================
# output info
#==============================================================
def output_info (cutoff, cluster, cluster_id, pfams, pfam_ann, assign_flag, outfile):
	details = {}
	details_rep = {}
	pers = {}
	cluster_flt = {}
	for myclust in cluster.keys():
		mytotal = len(cluster[myclust].keys())
		number = {}
		for myid in cluster[myclust].keys():
			if myid in pfams:
				for item in pfams[myid].keys():
					if not item in number:
						number[item] = {}
					number[item][myid] = ""
					#print(myclust + "\t" + item + "\t" + str(number[item].keys()))
        # foreach member
		if len(number.keys()) == 0:
			continue
		
		cluster_flt[myclust] = ""
		for item in number.keys():
			mynum = len(number[item].keys())
			myper = round(float(mynum)/float(mytotal), 2)
			if not myper in pers:
				pers[myper] = {}
			if not mytotal in pers[myper]:
				pers[myper][mytotal] = {}
			pers[myper][mytotal][myclust] = ""
			if myper >= float(cutoff):
				if not myclust in details:
					details[myclust] = {}
				if not item in details[myclust]:
					details[myclust][item] = myper
			if myclust in pfams:
				if not myclust in details_rep:
					details_rep[myclust] = {}
				if item in pfams[myclust]:
					details_rep[myclust][item] = myper
    # foreach cluster

	outfile1 = re.sub(".tsv", ".split.tsv", outfile)
	open_file = open(outfile, "w")
	open_file1 = open(outfile1, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tdescription\tconsistency"
	open_file.write(title + "\n")
	open_file1.write(title + "\n")
	for myclust in sorted(cluster_flt.keys()):
		clust_id = cluster_id[myclust]
		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in details_rep:
				ann_info = ""
				score_info = ""
				for item in sorted(details_rep[myclust].keys()):
					myann = "NA"
					if item in pfam_ann:
						ann_info = ann_info + pfam_ann[item] + ";"
						myann = pfam_ann[item]
					else:
						ann_info = ann_info + "NA;"
					score_info = score_info + str(details_rep[myclust][item]) + ";"
					open_file1.write(clust_id + "\tPfam_PfamDomain\t" + item + "\t" + str(details_rep[myclust][item]) + "\n")
				ann_info = re.sub(";$", "", ann_info)
				score_info = re.sub(";$", "", score_info)
				mypfam = ";".join(sorted(details_rep[myclust].keys()))
				open_file.write(clust_id + "\tPfam_PfamDomain\t" + mypfam + "\t" + ann_info + "\t" + score_info + "\n")
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in details:
				ann_info = ""
				score_info = ""
				for item in sorted(details[myclust].keys()):
					myann = "NA"
					if item in pfam_ann:
						ann_info = ann_info + pfam_ann[item] + ";"
						myann = pfam_ann[item]
					else:
						ann_info = ann_info + "NA;"
					score_info = score_info + str(details[myclust][item]) + ";"
					open_file1.write(clust_id + "\tPfam_PfamDomain\t" + item + "\t" + myann + "\t" + str(details[myclust][item]) + "\n")
				ann_info = re.sub(";$", "", ann_info)
				score_info = re.sub(";$", "", score_info)
				mypfam = ";".join(sorted(details[myclust].keys()))
				open_file.write(clust_id + "\tPfam_PfamDomain\t" + mypfam + "\t" + ann_info + "\t" + score_info + "\n")
		# different assignment method
	# foreach cluster
	open_file.close()
	
	outfile2 = re.sub(".tsv", ".all.spectrum.tsv", outfile)
	open_file2 = open(outfile2, "w")
	title = "percentage\tmember_num\tnumber"
	open_file2.write(title + "\n")
	for myper in sorted(pers.keys(), key=float):
		for mem_num in sorted(pers[myper], key=int):
			mynum = len(pers[myper][mem_num].keys())
			open_file2.write(str(myper) + "\t" + str(mem_num) + "\t" + str(mynum) + "\n")
    # foreach percentage
	open_file2.close()
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start interproscan_pfam_protein_family.py -p " + values.path + " ####\n")
	
	### collect sample info ###
	sys.stderr.write("Get info ......starting\n")
	cluster, cluster_id, cluster_mem = collect_cluster_info (config.protein_family)

	### collect annotation info and do stat ###
	sys.stderr.write("Get pfam info ......starting\n")
	pfam_ann = collect_pfam_ann (config.pfam_database)
	pfams = collect_pfam_info (cluster_mem, values.extension, values.path, values.output)
	sys.stderr.write("Get pfam info ......done\n")
	
	### Output sample info
	sys.stderr.write("\nOutput Pfam summary info ......starting\n")
	output_info (config.tshld_consistency, cluster, cluster_id, pfams, pfam_ann, values.method, values.output)
	sys.stderr.write("Output Pfam summary info ......done\n")

	sys.stderr.write("### Finish interproscan_pfam_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

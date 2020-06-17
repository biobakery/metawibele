#!/usr/bin/env python

"""
MetaWIBELE: interproscan_protein_family module
Summary the results of annotation from InterProScan

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
Summary the results of annotation from InterProScan
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
	                    default="interproscan.txt")
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
	cluster_id = {}
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
# collect annotation info
#==============================================================
def collect_ann_info (cluster_mem, extension, ann_path, types, outfile):	# list.txt
	anns = {}
	anns_info = {}
	filelist = utilities.find_files(ann_path, extension, None)
	for myfile in filelist:
		for suffix in types:
			myfile1 = re.sub(extension, suffix, myfile)
			#myfile = ann_path + "/" + samplelist + "/" + samplelist + "." + suffix # split1.interpro.SUPERFAMILY.tsv
			if not os.path.isfile(myfile1):
				print ("File not exist!\t" + myfile1)
				continue
			open_file = open(myfile1, "r")
			mym = re.search("interpro.([\S]+).tsv", suffix)
			mytype = "InterProScan_" + mym.group(1)
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
				acc = info[1]
				if not myid in anns:
					anns[myid] = {}
				if not mytype in anns[myid]:
					anns[myid][mytype] = {}
				if not mytype in anns_info:
					anns_info[mytype] = {}
				if info[2] == "":
					info[2] = "NA"
				if info[3] == "":
					info[3] = "NA"
				anns[myid][mytype][acc] = info[2] + "\t" + info[3]
				anns_info[mytype][acc] = info[2] + "\t" + info[3]
			# foreach line
			open_file.close()
		# foreach type of annotation
	#foreach samplelist

	# output details
	outfile1 = re.sub("_proteinfamilies.", "_proteinfamilies.ORF.", outfile)
	open_out = open(outfile1, "w")
	open_out.write(utilities.PROTEIN_ID + "\ttype\tdetail\tdescription\tInterPro_accession\n")
	for myid in sorted(anns.keys()):
		for mytype in sorted(anns[myid].keys()):
			myinfo1 = ""
			myinfo2 = ""
			myinfo3 = ""
			for myacc in sorted(anns[myid][mytype].keys()):
				myinfo1 = myinfo1 + myacc + ";"
				tmp = anns[myid][mytype][myacc].split("\t")
				myinfo2 = myinfo2 + tmp[0] + ";"
				myinfo3 = myinfo3 + tmp[1] + ";"
			# foreach annotated accession
			myinfo1 = re.sub(";$", "", myinfo1)
			myinfo2 = re.sub(";$", "", myinfo2)
			myinfo3 = re.sub(";$", "", myinfo3)
			if myinfo1 == "":
				continue
			open_out.write(myid + "\t" + mytype + "\t" + myinfo1 + "\t" + myinfo2 + "\t" + myinfo3 + "\n")
		# foreach type
	# foreach seqID
	open_out.close()
	return anns, anns_info
# function collect_ann_info


#==============================================================
# output info
#==============================================================
def output_info (cutoff, cluster, cluster_id, anns, anns_info, assign_flag, outfile):
	details = {}
	details_rep = {}
	cluster_flt = {}
	pers = {}
	for myclust in cluster.keys():
		mytotal = len(cluster[myclust].keys())
		number = {}
		for myid in cluster[myclust].keys():
			if myid in anns:
				for mytype in anns[myid].keys():
					if not mytype in number:
						number[mytype] = {}
					for item in anns[myid][mytype].keys():
						if not item in number[mytype]:
							number[mytype][item] = {}
						number[mytype][item][myid] = ""
					# debug
					#print(mytype + "\t" + str(number[mytype].keys()))
        # foreach member
		for mytype in number.keys():
			cluster_flt[myclust] = ""
			for item in number[mytype].keys():
				mynum = len(number[mytype][item].keys())
				myper = round(float(mynum)/float(mytotal), 2)
				# debug
				#print(mytype + "\t" + myclust + "\t" + item)
				if not mytype in pers:
					pers[mytype] = {}
				if not myper in pers[mytype]:
					pers[mytype][myper] = {}
				if not mytotal in pers[mytype][myper]:
					pers[mytype][myper][mytotal] =  {}
				pers[mytype][myper][mytotal][myclust] = ""
				if myper >= float(cutoff):
					if not myclust in details:
						details[myclust] = {}
					if not mytype in details[myclust]:
						details[myclust][mytype] = {}
					if not item in details[myclust][mytype]:
						details[myclust][mytype][item] = myper
				if myclust in anns:
					if not myclust in details_rep:
						details_rep[myclust] = {}
					if mytype in anns[myclust]:
						if not mytype in details_rep[myclust]:
							details_rep[myclust][mytype] = {}
						if item in anns[myclust][mytype]:
							details_rep[myclust][mytype][item] = myper
			# foreach item
		# foreach type
	# foreach cluster	

	# output info
	outfile1 = re.sub(".tsv", ".split.tsv", outfile)
	open_file = open(outfile, "w")
	open_file1 = open(outfile1, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tdescription\tInterPro_accession\tconsistency"
	open_file.write(title + "\n")
	open_file1.write(title + "\n")
	for myclust in sorted(cluster_flt.keys()):
		clust_id = cluster_id[myclust]
		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in details_rep:
				for mytype in sorted(details_rep[myclust].keys()):
					myacc = ""
					myann = ""
					myinter = ""
					myscore = ""
					for item in sorted(details_rep[myclust][mytype].keys()):
						if myacc == "":
							myacc = item
						else:
							myacc = myacc + ";" + item
						if mytype in anns_info:
							if item in anns_info[mytype]:
								ann_tmp, inter_tmp = anns_info[mytype][item].split("\t")
							else:
								ann_tmp = "NA"
								inter_tmp = "NA"
						if myann == "":
							myann = ann_tmp
						else:
							myann = myann + ";" + ann_tmp
						if myinter == "":
							myinter = inter_tmp
						else:
							myinter = myinter + ";" + inter_tmp
						if myscore == "":
							myscore = str(details_rep[myclust][mytype][item])
						else:
							myscore = myscore + ";" + str(details_rep[myclust][mytype][item])
						open_file1.write(clust_id + "\t" + mytype + "\t" + item + "\t" + ann_tmp + "\t" + inter_tmp + "\t" + str(details_rep[myclust][mytype][item]) + "\n")
					open_file.write(clust_id + "\t" + mytype + "\t" + myacc + "\t" + myann + "\t" + myinter + "\t" + myscore + "\n")
				# foreach type
			# if details_rep
		# if centroid
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in details:
				for mytype in sorted(details[myclust].keys()):
					myacc = ""
					myann = ""
					myinter = ""
					myscore = ""
					for item in sorted(details[myclust][mytype].keys()):
						if myacc == "":
							myacc = item
						else:
							myacc = myacc + ";" + item
						if mytype in anns_info:
							if item in anns_info[mytype]:
								ann_tmp, inter_tmp = anns_info[mytype][item].split("\t")
							else:
								ann_tmp = "NA"
								inter_tmp = "NA"
						if myann == "":
							myann = ann_tmp
						else:
							myann = myann + ";" + ann_tmp
						if myinter == "":
							myinter = inter_tmp
						else:
							myinter = myinter + ";" + inter_tmp
						if myscore == "":
							myscore = str(details[myclust][mytype][item])
						else:
							myscore = myscore + ";" + str(details[myclust][mytype][item])
						open_file1.write(clust_id + "\t" + mytype + "\t" + item + "\t" + ann_tmp + "\t" + inter_tmp + "\t" + str(details[myclust][mytype][item]) + "\n")
					open_file.write(clust_id + "\t" + mytype + "\t" + myacc + "\t" + myann + "\t" + myinter + "\t" + myscore + "\n")
				# foreach type
			# if details
		# if consistency
	# foreach cluster
	open_file.close()

	for mytype in sorted(pers.keys()):
		outfile2 = re.sub("detail.tsv", mytype + ".all.spectrum.tsv", outfile)
		open_file2 = open(outfile2, "w")
		title = "percentage\tmember_num\tnumber"
		open_file2.write(title + "\n")
		for myper in sorted(pers[mytype].keys(), key=float):
			for mem_num in sorted(pers[mytype][myper], key=int):
				mynum = len(pers[mytype][myper][mem_num].keys())
				open_file2.write(str(myper) + "\t" + str(mem_num) + "\t" + str(mynum) + "\n")
		open_file2.close()
	# foreach mytype
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start interproscan_protein_family.py -p " + values.path + " ####\n")
	
	### collect sample info ###
	sys.stderr.write("Get info ......starting\n")
	cluster, cluster_id, cluster_mem = collect_cluster_info (config.protein_family)

	### collect annotation info and do stat ###
	sys.stderr.write("Get annotation info ......starting\n")
	anns, anns_info = collect_ann_info (cluster_mem, values.extension, values.path, config.interproscan_type, values.output)
	sys.stderr.write("Get annotation info ......done\n")
	
	### Output sample info
	sys.stderr.write("\nOutput annotation summary info ......starting\n")
	output_info (config.tshld_consistency, cluster, cluster_id, anns, anns_info, values.method, values.output)
	sys.stderr.write("Output annotation summary info ......done\n")

	sys.stderr.write("### Finish interproscan_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python

"""
MetaWIBELE: ddi_DOMINE_protein_family module
Summary the results of DDI for families

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
	                    default="interpro.DDI.tsv")
	parser.add_argument('-p', "--path",
	                    help='input the path of annotation file',
	                    required=True)
	parser.add_argument('-a', "--method",
	                    help='specify how to assign annotations to families',
	                    choices=["centroid", "consistency"],
	                    required=True,
	                    default="consistency")
	parser.add_argument('-l', "--label",
	                    help='spefify the label name for DDI',
	                    default="DOMINE_interaction")
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
# collect DDI info
#==============================================================
def collect_DDI_info (cluster_mem, extension, ann_path, level, label, outfile):	# list.txt
	DDIs = {}
	anns = {}
	titles = {}
	filelist = utilities.find_files(ann_path, extension, None)
	for myfile in filelist:
		#myfile = ann_path + "/" + samplelist + "/" + samplelist + "." + suffix
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
			continue
		open_file = open(myfile, "r")
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
			mytype = info[titles["Type"]]
			mylevel = info[titles["Interaction"]]
			mypfam = info[titles["Pfam1_ID"]] + ":" + info[titles["Pfam2_ID"]]
			myann = info[titles["Pfam1_ann"]] + ":" + info[titles["Pfam2_ann"]]
			if not myid in DDIs:
				DDIs[myid] = {}
			DDIs[myid][mylevel + "\t" + mypfam] = ""
			anns[mypfam] = myann
		# foreach line
		open_file.close()
	# foreach samplelist

	# output details
	outfile1 = re.sub("_proteinfamilies.", "_proteinfamilies.ORF.", outfile)
	outfile2 = re.sub(".tsv", ".detail.tsv", outfile1)
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile2, "w")
	open_out2.write(utilities.PROTEIN_ID + "\ttype\tdetail\tannotation\tinteraction\n")
	open_out1.write(utilities.PROTEIN_ID + "\tType\tInteraction\tPfam1_ID\tPfam2_ID\tPfam1_ann\tPfam2_ann\n")
	for myid in sorted(DDIs.keys()):
		mypfam = ""
		myann = ""
		mylevel = ""
		for item in sorted(DDIs[myid].keys()):
			tmp = item.split("\t")
			#myt = "DOMINE_interaction"
			myt = label
			myl = tmp[0]
			pfam1, pfam2 = tmp[1].split(":")
			ann1 = "NA:NA"
			if tmp[1] in anns:
				ann1 = anns[tmp[1]]
			ann1 = re.sub(":", "\t", ann1)
			open_out1.write(myid + "\t" + myt + "\t" + myl + "\t" + pfam1 + "\t" + pfam2 + "\t" + ann1 + "\n")
			if level != "no":
				if tmp[0] != "NA":
					if tmp[0] != level:
						continue
			mypfam = mypfam + tmp[1] + ";"
			mylevel = mylevel + tmp[0] + ";"
			if tmp[1] in anns:
				myann = myann + anns[tmp[1]] + ";"
			else:
				myann = myann + "NA;"
		# foreach DDI
		mypfam = re.sub(";$", "", mypfam)
		myann = re.sub(";$", "", myann)
		mylevel = re.sub(";$", "", mylevel)
		if mypfam == "":
			continue
		#open_out2.write(myid + "\tDOMINE_interaction\t" + mypfam + "\t" + myann  + "\t" + mylevel + "\n")
		open_out2.write(myid + "\t" + label + "\t" + mypfam + "\t" + myann  + "\t" + mylevel + "\n")
	# foreach seqID
	open_out1.close()
	open_out2.close()

	return DDIs, anns
# function collect_DDI_info


#==============================================================
# output info
#==============================================================
def output_info (cutoff, cluster, cluster_id, DDIs, anns, level, assign_flag, label, outfile):
	details = {}
	details_rep = {}
	pers = {}
	for myclust in cluster.keys():
		mytotal = len(cluster[myclust].keys())
		number = {}
		for myid in cluster[myclust].keys():
			if myid in DDIs:
				# debug
				#print("Annotated member\t" + myclust + "\t" + myid)
				for item in DDIs[myid].keys():
					if not item in number:
						number[item] = {}
					number[item][myid] = ""
					#print(myclust + "\t" + item + "\t" + str(number[item].keys()))
        # foreach member
		for item in number.keys():
			mynum = len(number[item].keys())
			myper = round(float(mynum)/float(mytotal), 2)
			# debug
			#print(myclust + "\t" + item + "\t" + str(myper))
			mylevel, mypfam = item.split("\t")
			if level != "no":
				if mylevel != "NA":
					if mylevel == level:
						if not myper in pers:
							pers[myper] = {}
						if not mytotal in pers[myper]:
							pers[myper][mytotal] = {}
						pers[myper][mytotal][myclust] = ""
			else:
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
			if myclust in DDIs:
				if item in DDIs[myclust]:
					if not myclust in details_rep:
						details_rep[myclust] = {}
					details_rep[myclust][item] = myper
    # foreach cluster

	# get refined details
	details_flt = {}
	outs_detail = {}
	for myclust in sorted(cluster.keys()):
		clust_id = cluster_id[myclust]
		if assign_flag == "centroid":  # assign annotation based on representative info
			if myclust in details_rep:
				if not myclust in outs_detial:
					outs_detail[myclust] = {}
				for item in details_rep[myclust].keys():
					outs_detail[myclust][item] = details_rep[myclust][item]
					mylevel, mypfam = item.split("\t")
					if level != "no":
						if mylevel != "NA":
							if mylevel != level:
								continue
					if not myclust in details_flt:
						details_flt[myclust] = {}
					details_flt[myclust][item] = details[myclust][item]
		if assign_flag == "consistency":  # assign annotation baed on consistency info
			if myclust in details:
				if not myclust in outs_detail:
					outs_detail[myclust] = {}
				for item in details[myclust].keys():
					outs_detail[myclust][item] = details[myclust][item]
					mylevel, mypfam = item.split("\t")
					if level != "no":
						if mylevel != "NA":
							if mylevel != level:
								continue
					if not myclust in details_flt:
						details_flt[myclust] = {}
					details_flt[myclust][item] = details[myclust][item]
	# foreach cluster 

	# get number
	number = {}
	total_num = 0 
	for myclust in cluster.keys():
		total_num = total_num + 1 
		if myclust in outs_detail:
			for item in outs_detail[myclust].keys():
				mylevel, myitem = item.split("\t")
				if not mylevel in number:
					number[mylevel] = {}
				number[mylevel][myclust] = ""
            # foreach level
        # if interact
		else:
			if not "No_interaction" in number:
				number["No_interaction"] = {}
			number["No_interaction"][myclust] = ""
    # foreach cluster

	# output info
	open_file = open(outfile, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\tType\tInteraction\tPfam1_ID\tPfam2_ID\tPfam1_ann\tPfam2_ann"
	open_file.write(title + "\n")
	for myclust in sorted(outs_detail.keys()):
		clust_id = cluster_id[myclust]
		for item in sorted(outs_detail[myclust].keys()):
			mylevel, mypfam = item.split("\t")
			myann = "NA:NA"
			if mypfam in anns:
				myann = anns[mypfam]
			mypfam = re.sub(":", "\t", mypfam)
			myann = re.sub(":", "\t", myann)
			#open_file.write(clust_id + "\tDOMINE_interaction\t" + mylevel + "\t" + mypfam  + "\t" + myann + "\n")
			open_file.write(clust_id + "\t" + label + "\t" + mylevel + "\t" + mypfam  + "\t" + myann + "\n")
	open_file.close()

	outfile2 = re.sub(".tsv", ".detail.tsv", outfile)
	outfile3 = re.sub(".tsv", ".detail.split.tsv", outfile)
	open_file = open(outfile2, "w")
	open_file1 = open(outfile3, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\ttype\tdetail\tdescription\tconsistency"
	open_file.write(title + "\n")
	open_file1.write(title + "\n")
	for myclust in sorted(details_flt.keys()):
		clust_id = cluster_id[myclust]
		ddi_info = ""
		ann_info = ""
		score_info = ""
		level_info = ""
		for item in sorted(details_flt[myclust].keys()):
			mylevel, myitem = item.split("\t")
			myann = "NA"
			if myitem in anns:
				myann = anns[myitem]
			ann_info = ann_info + myann + ";"
			ddi_info = ddi_info + myitem + ";"
			level_info = level_info + mylevel + ";"
			score_info = score_info + str(details_flt[myclust][item]) + ";"
			#open_file1.write(clust_id + "\tDOMINE_interaction\t" + myitem + "\t" + myann + "\t" + str(details_flt[myclust][item]) + "\n")
			open_file1.write(clust_id + "\t" + label + "\t" + myitem + "\t" + myann + "\t" + str(details_flt[myclust][item]) + "\n")
		ddi_info = re.sub(";$", "", ddi_info)
		ann_info = re.sub(";$", "", ann_info)
		score_info = re.sub(";$", "", score_info)
		level_info = re.sub(";$", "", level_info)
		#open_file.write(clust_id + "\tDOMINE_interaction\t" + ddi_info + "\t" + ann_info + "\t" + score_info + "\n")
		open_file.write(clust_id + "\t" + label + "\t" + ddi_info + "\t" + ann_info + "\t" + score_info + "\n")
	open_file.close()
	open_file1.close()
	
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

	sys.stderr.write("### Start ddi_DOMINE_protein_family.py -p " + values.path + " ####\n")
	
	### collect sample info ###
	sys.stderr.write("Get info ......starting\n")
	cluster, cluster_id, cluster_mem = collect_cluster_info (config.protein_family)

	### collect annotation info and do stat ###
	sys.stderr.write("Get pfam info ......starting\n")
	DDIs, anns = collect_DDI_info (cluster_mem, values.extension, values.path, "HC", values.label, values.output)
	sys.stderr.write("Get pfam info ......done\n")
	
	### Output sample info
	sys.stderr.write("\nOutput Pfam summary info ......starting\n")
	output_info (config.tshld_consistency, cluster, cluster_id, DDIs, anns, "HC", values.method, values.label, values.output)
	sys.stderr.write("Output Pfam summary info ......done\n")

	sys.stderr.write("### Finish ddi_DOMINE_protein_family.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

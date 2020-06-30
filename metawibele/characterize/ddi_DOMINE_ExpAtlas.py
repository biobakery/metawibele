#!/usr/bin/env python

"""
MetaWIBELE: ddi_DOMINE_ExpAtlas module
Based on Expression Atlas DB, check whether human transcripts expressed in gut for each ORF

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
Based on Expression Atlas DB, check whether human transcripts expressed in gut
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--ddi-annotation",
	                    help='DDIs annotation file',
	                    required=True)
	parser.add_argument('-o', "--output",
	                    help='output summary ann file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect info
#==============================================================
def collect_expression_info (exp_list):	#  expression.list
	expression = {}
	for myfile in exp_list:
		myfile = myfile.strip()
		if not len(myfile):
			continue
		if re.search("^#", myfile):
			continue
		for line in utils.gzip_bzip2_biom_open_readlines (myfile): 
			line = line.strip()
			if not len(line):
				continue
			if re.search("^#", line) or re.search("^Gene ID", line):
				continue
			info = line.split("\t")
			#flag = 0
			#myindex = 2
			#while myindex < len(info):
			#	if info[myindex] == "":
			#		continue
			#	if float(info[myindex]) >= float(cutoff):
			#		flag = 1
			#	myindex = myindex + 1
			#if flag == 1:
			mygene = info[1]
			expression[mygene] = ""
		# foreach line
	# foreach dataset
	
	return expression
# function collect_expression_info


def collect_pfam_info (map_file):	# uniprot_human_pfam.tsv 
	pfams = {}
	titles = {}
	for line in utils.gzip_bzip2_biom_open_readlines (map_file):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		if re.search("^Pfam\t", line):
			info = line.split("\t")
			myindex = 0
			while myindex < len(info):
				titles[info[myindex]] = myindex
				myindex = myindex + 1
			continue
		info = line.split("\t")
		if len(info) != len(titles.keys()):
			continue
		pfam = info[titles["Pfam"]]
		gene = info[titles["Gene_names"]]
		taxaID = info[titles["NCBI_TaxID"]]
		if taxaID != "9606":	# filter out for human pfam
			continue
		pfams[pfam] = gene
	# foreach line
	
	return pfams
# function collect_pfam_info


def collect_interaction_info (ann_file):	# summary_DOMINE_peptide.ann.tsv
	interact = {}
	titles = {}
	title = ""
	open_file = open(ann_file, "r")
	line = open_file.readline()
	line = line.strip()
	title = line
	info = line.split("\t")
	myindex = 0
	myflag = info[0]
	while myindex < len(info):
		titles[info[myindex]] = myindex
		myindex = myindex + 1
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid1 = info[0]
		pfam1 = info[titles["Pfam1_ID"]]
		pfam2 = info[titles["Pfam2_ID"]]
		myid2 = pfam1 + "\t" + pfam2
		if not myid1 in interact:
			interact[myid1] = {}
		interact[myid1][myid2] = line
	# foreach line
	open_file.close()
	return interact, title, myflag
# collect_interaction_info


#==============================================================
# assign Expression Atlas supporting info to DDIs 
#==============================================================
def assign_interaction (cutoff, pfams, expression, interact, title, id_flag, outfile):
	# check expression
	exp = {}
	for mypfam in pfams.keys():	# check human Pfam
		info = pfams[mypfam].split(";")
		gene_total = len(info)
		exp_num = 0
		for item in info:
			if item in expression:
				exp_num = exp_num + 1
		# foreach gene
		per = 1.0 * exp_num / gene_total
		if exp_num > 0:
			exp[mypfam] = pfams[mypfam] + "\t" + str(per)
	# foreach pfam

	# check interaction
	number = {}
	anns = {}
	pers_all = {}
	for pepid in sorted(interact.keys()):
		flag = 0
		for pfam_pair in sorted(interact[pepid].keys()):
			pfam1, pfam2 = pfam_pair.split("\t")	
			if pfam1 in pfams or pfam2 in pfams:	# have human Pfam signature
				flag = 1
				myinter = "ExpAtlas_human_unknown"
				score = 0
				if pfam1 in pfams and pfam1 in exp:	# human gene including this domain expressed in human hut
					gene_info, per = exp[pfam1].split("\t")
					if not pepid in pers_all:
						pers_all[pepid] = {}
					if not pfam_pair in pers_all[pepid]:
						pers_all[pepid][pfam_pair] = round(float(per), 2)
					if float(per) >= float(cutoff):
						flag = 2
						myinter = "ExpAtlas_human_support"
						score = per
				if flag == 1:
					if pfam2 in pfams and pfam2 in exp:	# human gene including this domain expressed in human hut
						gene_info, per = exp[pfam2].split("\t")
						if not pepid in pers_all:
							pers_all[pepid] = {}
						if not pfam_pair in pers_all[pepid]:
							pers_all[pepid][pfam_pair] = round(float(per), 2)
						if float(per) > float(cutoff):
							myinter = "ExpAtlas_human_support"
							score = per
							flag = 2
				if flag == 2:
					if not pepid in anns:
						anns[pepid] = {}
					if not pfam_pair in anns[pepid]:
						anns[pepid][pfam_pair] = interact[pepid][pfam_pair] + "\t" + myinter + "\t" + str(score)
				if not myinter in number:
					number[myinter] = {}
				number[myinter][pepid] = pfam_pair
		# foreach pair
		if flag == 0: # no human Pfam hits
			myinter = "No_human_pfam"
			if not myinter in number:
				number[myinter] = {}
			number[myinter][pepid] = pfam_pair
	# foreach cluster

	outs = {}
	details = {}
	outfile1 = re.sub(".detail.tsv" ,".tsv", outfile)
	open_out1 = open(outfile1, "w")
	open_out2 = open(outfile, "w")
	open_out1.write(title + "\t" + "ExpAtlas_type\tExpAtlas_score" + "\n")
	open_out2.write(id_flag + "\ttype\tdetail\tExpAtlas\n")
	for mypep in sorted(anns.keys()):
		for myinfo in sorted(anns[mypep].keys()):
			ddi = re.sub("\t", ":", myinfo)
			info = anns[mypep][myinfo].split("\t")
			myscore = info[-1]
			myann = info[1] + "\t" + info[2] + "\t" + info[3] + "\t" + info[4]
			if not myann in outs:
				outs[myann] = {}
			outs[myann][mypep] = ""
			if not mypep in details:
				details[mypep] = {}
			details[mypep][ddi] = myscore
			open_out1.write(anns[mypep][myinfo] + "\n")
		# foreach item
	# foreah pepid
	open_out1.close()
	for mypep in sorted(details.keys()):
		ddi_str = ""
		score_str = ""
		for myddi in sorted(details[mypep].keys()):
			ddi_str = ddi_str + myddi + ";"
			score_str = score_str + details[mypep][myddi] + ";"
		# foreach DDI
		ddi_str = re.sub(";$", "", ddi_str)
		score_str = re.sub(";$", "", score_str)
		open_out2.write(mypep + "\tExpAtlas_interaction\t" + ddi_str + "\t" + score_str + "\n")
	# foreach cluster
	open_out2.close()

	outfile3 = re.sub(".tsv", ".info.tsv", outfile1)
	open_out3 = open(outfile3, "w")
	open_out3.write("Pfam1_ID\tPfam1_ann\tPfam2_ID\tPfam2_ann\tnumber\n")
	orders = {}
	for myid in outs.keys():
		mynum = len(outs[myid].keys())
		if not mynum in orders:
			orders[mynum] = {}
		orders[mynum][myid] = ""
	# foreach pair
	for mynum in sorted(orders.keys(), key=int, reverse=True):
		for myid in sorted(orders[mynum].keys()):
			open_out3.write(myid + "\t" + str(mynum) + "\n")
		# foreach id
	# foreach number
	open_out3.close()

# assign_annotation


#==============================================================
###########  Main processing ############
#=============================================================
def main():

	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start ddi_DOMINE_ExpAtlas.py -i " + values.ddi_annotation + " ####\n")
	

	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	pfams = collect_pfam_info (config.human_pfam_database)
	exp = collect_expression_info (config.Expression_Atlas_database)
	interact, title, id_flag = collect_interaction_info (values.ddi_annotation)
	sys.stderr.write("Get info ......done\n")
	
	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign ExpAtlas expression to peptide families ......starting\n")
	cutoff = 0
	assign_interaction (cutoff, pfams, exp, interact, title, id_flag, values.output)
	sys.stderr.write("\nAssign ExpAtlas expression to peptide families ......done\n")

	sys.stderr.write("### Finish ddi_DOMINE_ExpAtlas.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

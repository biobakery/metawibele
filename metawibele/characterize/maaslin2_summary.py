#!/usr/bin/env python

"""
MetaWIBELE: maaslin2_summary module
Summary association summary for protein families

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
import math

try:
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Summary association summary for protein families
"""

def get_args():

	parser=argparse.ArgumentParser()
	parser.add_argument('-a', "--abundance",
	                    help='input stat abundance info file',
	                    required=True)
	parser.add_argument('-b', "--fold-change",
	                    help='input fold change info file',
	                    required=True)
	parser.add_argument('-c', "--prevalence",
	                    help='input prevalence info file',
	                    required=True)
	parser.add_argument('-p', "--prevalence-cutoff",
	                    help='cutoff for prevalence filtering (prevalence >= cutoff), e.g. no|0.10',
	                    default="0.10")
	parser.add_argument('-q', "--qvalue-cutoff",
	                    help='cutoff for FDR adjust pvalue filtering (qvalue < cutoff), e.g. no|0.05',
	                    default="0.05")
	parser.add_argument('-e', "--effect-size",
	                    help='specify the item name indicating effect size',
						required=True)
	parser.add_argument('-o', "--output",
	                    help='output DA summary file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect differential abundance stat info
#==============================================================
def collect_DA_stat_info (stat_file):
	stat = {}
	# collect stat info
	sys.stderr.write("Get DA stat info ......\n")
	open_file = open(stat_file, "r")
	titles = {}
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_FAMILY_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myclust = info[titles[utilities.PROTEIN_FAMILY_ID]]
		mytype = info[titles["cmp_type"]]
		coef = info[titles["coef"]]
		stderr = info[titles["stderr"]]
		pvalue = info[titles["pval"]]
		qvalue = info[titles["qval"]]
		myid = myclust + "\t" + mytype
		stat[myid] = coef + "\t" + stderr + "\t" + pvalue + "\t" + qvalue
	# foreach line
	open_file.close()
	sys.stderr.write("Get DA stat info ......done\n")
	return stat
# collect_DA_stat_info


#==============================================================
# collect differential abundance prevalence info
#==============================================================
def collect_DA_prevalence_info (pre_file):
	prevalence = {}
	# collect prevalence info
	sys.stderr.write("Get DA prevalence info ......\n")
	open_file = open(pre_file, "r")
	titles = {}
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_FAMILY_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myclust = info[titles[utilities.PROTEIN_FAMILY_ID]]
		mytype = info[titles["cmp_type"]]
		pre = info[titles["prevalence"]]
		pre_yes = info[titles["prevalence_case"]]
		pre_no = info[titles["prevalence_control"]]
		myid = myclust + "\t" + mytype
		prevalence[myid] = pre + "\t" + pre_yes + "\t" + pre_no
	# foreach line
	open_file.close()
	sys.stderr.write("Get DA prevalence info ......done\n")
	return prevalence
# collect_DA_prevalence_info


#==============================================================
# collect fold change info 
#==============================================================
def collect_fold_change_info (fold_file, stat, prevalence, effect_size):
	# collect fold info
	folds = {}
	titles = {}
	sys.stderr.write("Get fold change info ......\n")
	open_file = open(fold_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_FAMILY_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myclust = info[titles[utilities.PROTEIN_FAMILY_ID]]
		mytype = info[titles["cmp_type"]]
		fold = info[titles["log(FC)"]]
		var = info[titles["mean(log)"]]
		effect = info[titles["Cohen's_d"]]
		myabun = info[titles["mean_abundance"]]
		myabun_case = info[titles["mean_abundance_case"]]
		myabun_control = info[titles["mean_abundance_control"]]
		myabun2 = info[titles["mean_prevalent_abundance"]]
		myabun2_case = info[titles["mean_prevalent_abundance_case"]]
		myabun2_control = info[titles["mean_prevalent_abundance_control"]]
		myid = myclust + "\t" + mytype
		mypvalue = "NA\tNA\tNA\tNA"
		mypre = "NA\tNA\tNA"
		if myid in stat:
			mypvalue = str(stat[myid])
		if myid in prevalence:
			mypre = prevalence[myid]
		if mypvalue != "NA\tNA\tNA\tNA":
			mycoef, stderr, myp, myq = mypvalue.split("\t")
			if myq == "NA" or myq == "NaN" or myq == "nan":
				continue
			mydis = abs(float(mycoef))
			if effect_size == "Cohen's_d":
				mydis = abs(float(effect))
			if effect_size == "log(FC)":
				mydis = abs(float(fold))
			if effect_size == "mean(log)":
				mydis = abs(float(var))
			if not mydis in folds:
				folds[mydis] = {}
			if not myq in folds[mydis]:
				folds[mydis][myq] = {}
			folds[mydis][myq][myid] = mypre + "\t" + mypvalue + "\t" + fold + "\t" + effect + "\t" + var + "\t" + str(myabun) + "\t" + str(myabun_case) + "\t" + str(myabun_control) + "\t" + str(myabun2) + "\t" + str(myabun2_case) + "\t" + str(myabun2_control)
	# foreach line
	open_file.close()
	sys.stderr.write("Get fold change info ......done\n")
	return folds
# collect_fold_change_info


#==============================================================
# summary differential abundance info
#==============================================================
def summary_info (folds, p_cutoff, q_value_cutoff, outfile):
	diff = {}
	types = {}
	diff_all = {}
	types_all = {}
	
	#### combine info and output info ####
	sys.stderr.write("Combine info and output info ......\n")
	outfile3 = re.sub(".tsv", ".all.tsv", outfile)
	pres = {}
	pres_all = {}
	pres_yes = {}
	pres_yes_all = {}
	pres_no = {}
	pres_no_all = {}
	qvalues = {}
	qvalues_all = {}
	abun = {}
	abun_all = {}
	open_file = open(outfile, "w")
	open_file3 = open(outfile3, "w")
	title = utilities.PROTEIN_FAMILY_ID + "\tcmp_type\tprevalence\tprevalence_case\tprevalence_control\tcoef\tstderr\tpvalue\tqvalue\tlog(FC)\tCohen's_d\tmean(log)\tmean_abundance\tmean_abundance_case\tmean_abundance_control\tmean_prevalent_abundance\tmean_prevalent_abundance_case\tmean_prevalent_abundance_control"
	open_file.write(title + "\n")
	open_file3.write(title + "\n")

	for mydis in sorted(folds.keys(), key=float, reverse=True):	# sort by effect size
		for myq in sorted(folds[mydis].keys(), key=float):	# sort by q-value
			for myid in sorted(folds[mydis][myq].keys()):
				open_file3.write(myid + "\t" + folds[mydis][myq][myid] + "\n")
				mypre, mypre_yes, mypre_no, mycoef, stderr, pvalue, qvalue, myfold, myeffect, myvar, myabun, myabun_case, myabun_control, myabun2, myabun2_case, myabun2_control = folds[mydis][myq][myid].split("\t")
				myclust, mytype = myid.split("\t")
				if not mytype in pres_all:
					pres_all[mytype] = {}
				mypre_tmp = round(float(mypre), 3)
				if not mypre_tmp in pres_all[mytype]:
					pres_all[mytype][mypre_tmp] = {}
				pres_all[mytype][mypre_tmp][myclust] = ""
				if not mytype in pres_yes_all:
					pres_yes_all[mytype] = {}
				mypre_yes_tmp = round(float(mypre_yes), 3)
				if not mypre_yes_tmp in pres_yes_all[mytype]:
					pres_yes_all[mytype][mypre_yes_tmp] = {}
				pres_yes_all[mytype][mypre_yes_tmp][myclust] = ""
				if not mytype in pres_no_all:
					pres_no_all[mytype] = {}
				mypre_no_tmp = round(float(mypre_no), 3)
				if not mypre_no_tmp in pres_no_all[mytype]:
					pres_no_all[mytype][mypre_no_tmp] = {}
				pres_no_all[mytype][mypre_no_tmp][myclust] = ""
				if not mytype in qvalues_all:
					qvalues_all[mytype] = {}
				qvalue_tmp = round(float(qvalue), 3)
				if not qvalue_tmp in qvalues_all[mytype]:
					qvalues_all[mytype][qvalue_tmp] = {}
				qvalues_all[mytype][qvalue_tmp][myclust] = ""
				types_all[mytype] = ""
				if not mytype in abun_all:
					abun_all[mytype] = {}
				abun_tmp = round(float(myabun), 3)
				if not abun_tmp in abun_all[mytype]:
					abun_all[mytype][abun_tmp] = {}
				abun_all[mytype][abun_tmp][myclust] = ""
				if not myclust in diff_all:
					diff_all[myclust] = {}
				if not mytype in diff_all[myclust]:
					diff_all[myclust][mytype] = folds[mydis][myq][myid] 

				### filtering
				flag = 1
				if p_cutoff != "no" and mypre != "NA":
					if float(mypre) < float(p_cutoff):	# filtering based on prevalence
						print("Filter out based on prevalence\t" + myid + "\t" + str(mypre))
						flag = 0	
				if q_value_cutoff != "no" and qvalue != "NA":
					if float(qvalue) >= float(q_value_cutoff):	# filtering based on q-value
						print("Filter out based on q-value\t" + myid + "\t" + str(qvalue))
						flag = 0	
				if flag == 1: # not filtering out
					if not mytype in pres:
						pres[mytype] = {}
					mypre_tmp = round(float(mypre), 3)
					if not mypre_tmp in pres[mytype]:
						pres[mytype][mypre_tmp] = {}
					pres[mytype][mypre_tmp][myclust] = ""
					if not mytype in pres_yes:
						pres_yes[mytype] = {}
					mypre_yes_tmp = round(float(mypre_yes), 3)
					if not mypre_yes_tmp in pres_yes[mytype]:
						pres_yes[mytype][mypre_yes_tmp] = {}
					pres_yes[mytype][mypre_yes_tmp][myclust] = ""
					if not mytype in pres_no:
						pres_no[mytype] = {}
					mypre_no_tmp = round(float(mypre_no), 3)
					if not mypre_no_tmp in pres_no[mytype]:
						pres_no[mytype][mypre_no_tmp] = {}
					pres_no[mytype][mypre_no_tmp][myclust] = ""
					if not mytype in qvalues:
						qvalues[mytype] = {}
					qvalue_tmp = round(float(qvalue), 3)
					if not qvalue_tmp in qvalues[mytype]:
						qvalues[mytype][qvalue_tmp] = {}
					qvalues[mytype][qvalue_tmp][myclust] = ""
					if not mytype in abun:
						abun[mytype] = {}
					abun_tmp = round(float(myabun), 3)
					if not abun_tmp in abun[mytype]:
						abun[mytype][abun_tmp] = {}
					abun[mytype][abun_tmp][myclust] = ""
					if not myclust in diff:
						diff[myclust] = {}
					if not mytype in diff[myclust]:
						diff[myclust][mytype] = folds[mydis][myq][myid] 
					types[mytype] = ""
					open_file.write(myclust + "\t" + mytype + "\t" + folds[mydis][myq][myid] + "\n")
			# foreach cluster
		# foreach abundance
	# foreach fold
	open_file.close()
	open_file3.close()

	# output prevalence distribution
	outfile4 = re.sub(".tsv", ".prevalence.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("prevalence\ttype\tnumber\n")
	for mytype in sorted(pres.keys()):
		for mypre in sorted(pres[mytype].keys(), key=float):
			mynum = 0
			for tmp_pre in sorted(pres[mytype].keys()):
				if float(tmp_pre) >= float(mypre):
					mynum = mynum + len(pres[mytype][tmp_pre])
			# foreach prevalence
			open_file.write(str(mypre) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach prevalence
	# foreach type
	open_file.close()
	outfile4 = re.sub(".tsv", ".all.prevalence.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("prevalence\ttype\tnumber\n")
	for mytype in sorted(pres_all.keys()):
		for mypre in sorted(pres_all[mytype].keys(), key=float):
			mynum = 0
			for tmp_pre in sorted(pres_all[mytype].keys()):
				if float(tmp_pre) >= float(mypre):
					mynum = mynum + len(pres_all[mytype][tmp_pre])
			# foreach prevalence
			open_file.write(str(mypre) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach prevalence
	# foreach type
	open_file.close()
	outfile4 = re.sub(".tsv", ".prevalence.case.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("prevalence\ttype\tnumber\n")
	for mytype in sorted(pres_yes.keys()):
		for mypre in sorted(pres_yes[mytype].keys(), key=float):
			mynum = 0
			for tmp_pre in sorted(pres_yes[mytype].keys()):
				if float(tmp_pre) >= float(mypre):
					mynum = mynum + len(pres_yes[mytype][tmp_pre])
			# foreach prevalence
			open_file.write(str(mypre) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach prevalence
	# foreach type
	open_file.close()
	outfile4 = re.sub(".tsv", ".all.prevalence.case.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("prevalence\ttype\tnumber\n")
	for mytype in sorted(pres_yes_all.keys()):
		for mypre in sorted(pres_yes_all[mytype].keys(), key=float):
			mynum = 0
			for tmp_pre in sorted(pres_yes_all[mytype].keys()):
				if float(tmp_pre) >= float(mypre):
					mynum = mynum + len(pres_yes_all[mytype][tmp_pre])
			# foreach prevalence
			open_file.write(str(mypre) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach prevalence
	# foreach type
	open_file.close()
	outfile4 = re.sub(".tsv", ".prevalence.control.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("prevalence\ttype\tnumber\n")
	for mytype in sorted(pres_no.keys()):
		for mypre in sorted(pres_no[mytype].keys(), key=float):
			mynum = 0
			for tmp_pre in sorted(pres_no[mytype].keys()):
				if float(tmp_pre) >= float(mypre):
					mynum = mynum + len(pres_no[mytype][tmp_pre])
			# foreach prevalence
			open_file.write(str(mypre) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach prevalence
	# foreach type
	open_file.close()
	outfile4 = re.sub(".tsv", ".all.prevalence.control.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("prevalence\ttype\tnumber\n")
	for mytype in sorted(pres_no_all.keys()):
		for mypre in sorted(pres_no_all[mytype].keys(), key=float):
			mynum = 0
			for tmp_pre in sorted(pres_no_all[mytype].keys()):
				if float(tmp_pre) >= float(mypre):
					mynum = mynum + len(pres_no_all[mytype][tmp_pre])
			# foreach prevalence
			open_file.write(str(mypre) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach prevalence
	# foreach type
	open_file.close()
	

	# output qvalue
	outfile4 = re.sub(".tsv", ".qvalue.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("qvalue\ttype\tnumber\n")
	for mytype in sorted(qvalues.keys()):
		for myq in sorted(qvalues[mytype].keys(), key=float):
			mynum = 0
			for tmp_q in sorted(qvalues[mytype].keys()):
				if float(tmp_q) <= float(myq):
					mynum = mynum + len(qvalues[mytype][tmp_q])
			# foreach qvalue
			open_file.write(str(myq) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach qvalue
	# foreach type
	open_file.close()
	outfile4 = re.sub(".tsv", ".all.qvalue.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("qvalue\ttype\tnumber\n")
	for mytype in sorted(qvalues_all.keys()):
		for myq in sorted(qvalues_all[mytype].keys(), key=float):
			mynum = 0
			for tmp_q in sorted(qvalues_all[mytype].keys()):
				if float(tmp_q) <= float(myq):
					mynum = mynum + len(qvalues_all[mytype][tmp_q])
			# foreach qvalue
			open_file.write(str(myq) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach qvalue
	# foreach type
	open_file.close()
	
	# output abundance distribution
	outfile4 = re.sub(".tsv", ".abundance.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("abundance\ttype\tnumber\n")
	for mytype in sorted(abun.keys()):
		for myab in sorted(abun[mytype].keys(), key=float):
			mynum = 0
			for tmp_ab in sorted(abun[mytype].keys()):
				if float(tmp_ab) >= float(myab):
					mynum = mynum + len(abun[mytype][tmp_ab])
			# foreach abundance
			open_file.write(str(myab) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach abundance
	# foreach type
	open_file.close()
	outfile4 = re.sub(".tsv", ".all.abundance.plot.tsv", outfile)
	open_file = open(outfile4, "w")
	open_file.write("abundance\ttype\tnumber\n")
	for mytype in sorted(abun_all.keys()):
		for myab in sorted(abun_all[mytype].keys(), key=float):
			mynum = 0
			for tmp_ab in sorted(abun_all[mytype].keys()):
				if float(tmp_ab) >= float(myab):
					mynum = mynum + len(abun_all[mytype][tmp_ab])
			# foreach abundance
			open_file.write(str(myab) + "\t" + mytype + "\t" + str(mynum) + "\n")
		# foreach abundance 
	# foreach type
	open_file.close()

	# venn plot
	outfile2 = re.sub(".tsv", ".venn.tsv", outfile)
	open_file2 = open(outfile2, "w")
	mystr = ""
	for mytype in sorted(types.keys()):
		mystr = mystr + mytype + "\t"
	# foreach type
	mystr = mystr.strip()
	open_file2.write(mystr + "\n")
	for myclust in sorted(diff.keys()):
		mystr = ""
		for mytype in sorted(types.keys()):
			if mytype in diff[myclust]:
				mystr = mystr + myclust + "\t"
			else:
				mystr = mystr + "NA\t"
		# foreach type
		mystr = mystr.strip()
		open_file2.write(mystr + "\n")
	# foreach cluster
	open_file2.close()
	outfile2 = re.sub(".tsv", ".all.venn.tsv", outfile)
	open_file2 = open(outfile2, "w")
	mystr = ""
	for mytype in sorted(types_all.keys()):
		mystr = mystr + mytype + "\t"
	# foreach type
	mystr = mystr.strip()
	open_file2.write(mystr + "\n")
	for myclust in sorted(diff_all.keys()):
		mystr = ""
		for mytype in sorted(types_all.keys()):
			if mytype in diff_all[myclust]:
				mystr = mystr + myclust + "\t"
			else:
				mystr = mystr + "NA\t"
		# foreach type
		mystr = mystr.strip()
		open_file2.write(mystr + "\n")
	# foreach cluster
	open_file2.close()
# function collect_diff_abundance_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args()
	values.effect_size = re.sub("mean_log", "mean(log)", values.effect_size)
	values.effect_size = re.sub("log_FC", "log(FC)", values.effect_size)

	sys.stderr.write("### Start maaslin2_summary.py -a " + values.abundance + " ####\n")
	

	### collect stat abundance info ###
	sys.stderr.write("Get stat abundance and fold change info ......starting\n")
	stat = collect_DA_stat_info (values.abundance)
	prevalence = collect_DA_prevalence_info (values.prevalence)
	folds = collect_fold_change_info (values.fold_change, stat, prevalence, values.effect_size)
	sys.stderr.write("Get stat abundance and fold change info ......done\n")
	
	### Output sample info
	sys.stderr.write("\nOutput diff abundance summary info ......starting\n")
	summary_info (folds, values.prevalence_cutoff, values.qvalue_cutoff, values.output)
	sys.stderr.write("Output diff abundance summary info ......done\n")

	sys.stderr.write("### Finish maaslin2_summary.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python3

"""
MetaWIBELE: maaslin2_collection module
Collect statistics info from MaAsLin2 results

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
	from metawibele import config
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Collect statistics info from MaAsLin2 results
"""


def get_args():

	parser = argparse.ArgumentParser()
	parser.add_argument('-a', "--abundance-table",
	                    help='input abundance table, features as rows and samples as columns',
	                    required=True)
	parser.add_argument('-s', "--smoothed-abundance",
	                    help='input smoothed abundance table, features as rows and samples as columns',
	                    required=True)
	parser.add_argument('-i', "--maaslin2-results",
	                    help='input the statistics results',
	                    required=True)
	parser.add_argument('-m', "--metadata",
	                    help='input the metadata',
	                    default="none")
	parser.add_argument('-o', "--output",
	                    help='output DA stat file',
	                    required=True)
	values = parser.parse_args()
	return values

# get_args


#==============================================================
# collect statistical results
#==============================================================
def collect_DA_info (DA_file):
	stats = {}
	titles = {}
	open_file = open(DA_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 0
	while myindex < len(info):
		item = info[myindex]
		titles[item] = myindex
		myindex = myindex + 1
	sys.stderr.write("Collect DA info......starting\n")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		meta = info[titles["metadata"]]
		value = info[titles["value"]]
		feature = info[titles["feature"]]
		if re.search("Cluster_[\d]+\.", feature):
			mym = re.search("(Cluster_[\d]+)", feature)
			myid = mym.group(1)
			feature = re.sub(myid + "\.", myid + "|", feature)
			feature = re.sub("\.\.", "@ ", feature)
			feature = re.sub("\.", " ", feature)
			feature = re.sub("\@", ".", feature)
		coef = info[titles["coef"]]
		stderr = info[titles["stderr"]]
		total = info[titles["N"]]
		non_zero = info[titles["N.not.0"]]
		p_value = info[titles["pval"]]
		q_value = info[titles["qvalue"]]
		if q_value == "NA":
			continue
		if meta in config.contrast_status:
			contrasts = config.contrast_status[meta].split(",")
			ref_status = config.ref_status[meta].split(",")
			tmp = {}
			for item in ref_status:
				mym = re.search("([\S]+)_vs_([\S]+)", item)
				mycon = mym.group(1)
				myref = mym.group(2)
				tmp[mycon] = myref
			if value in contrasts:
				if not value in tmp:
					continue
				myref = tmp[value]
				cmp_type = value + "_vs_" + myref
				prevalence = 1.0 * float(non_zero) / float(total)
				if not feature in stats:
					stats[feature] = {}
				stats[feature][cmp_type] = str(prevalence) + "\t" + coef + "\t" + stderr + "\t" + p_value + "\t" + q_value
	# foreach line
	open_file.close()
	sys.stderr.write("Collect DA info......done\n")
	return stats
# collect_DA_info



#==============================================================
# collect metadata for each sample
#==============================================================
def collect_meta_info (meta_file):
	meta_case = {}
	meta_control = {}
	titles = {}
	open_file = open(meta_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	metas = config.ref_status.keys()
	sys.stderr.write("Collect metadata info......starting\n")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		sample = info[0]
		if re.search("^[\d]+$", sample):
			sample = config.study + "_" + sample
		mymeta = ""
		for item in metas:
			myindex = titles[item]
			myitem = info[myindex]
			myref = config.ref_status[item].split(",")
			myrefs = {}
			mycons = {}
			for x in myref:
				mym = re.search("([\S]+)_vs_([\S]+)", x)
				mycon_tmp = mym.group(1)
				myref_tmp = mym.group(2)
				if not myref_tmp in myrefs:
					myrefs[myref_tmp] = []
				myrefs[myref_tmp].append(x)
				if not mycon_tmp in mycons:
					mycons[mycon_tmp] = []
				mycons[mycon_tmp].append(x)
			mycontrast = config.contrast_status[item].split(",")
			if myitem in myrefs:  # ref sample
				for mycmp in myrefs[myitem]:
					if not mycmp in meta_control:
						meta_control[mycmp] = {}
					meta_control[mycmp][sample] = ""
			if myitem in mycons: # contrast sample
				for mycmp in mycons[myitem]:
					if not mycmp in meta_case:
						meta_case[mycmp] = {}
					meta_case[mycmp][sample] = ""
	# foreach line
	open_file.close()
	sys.stderr.write("Collect metadata info......done\n")
	return meta_case, meta_control
# collect_meta_info


# ==============================================================
# collect abundance info and fold change
# ============================ =================================
def collect_abundance_info (abundance_file, DA_cluster, meta_case, meta_control):
	samples = {}
	total_value = {}
	abundance = {}
	meta_case_new = {}
	meta_control_new = {}
	open_file = open(abundance_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 1
	sys.stderr.write("Collect abundance info......starting\n")
	while myindex < len(info):
		item = info[myindex]
		samples[myindex] = item
		myindex = myindex + 1
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		cluster = info[0]
		if not cluster in DA_cluster:
			continue
		cmp_types = DA_cluster[cluster].keys()
		myindex = 1
		while myindex < len(info):
			item = info[myindex]
			mys = samples[myindex]
			for mycmp in cmp_types:
				if mycmp in meta_case:
					if mys in meta_case[mycmp]:
						if not cluster in total_value:
							total_value[cluster] = {}
						if not mycmp in total_value[cluster]:
							total_value[cluster][mycmp] = {}
						if not "yes" in total_value[cluster][mycmp]:
							total_value[cluster][mycmp]["yes"] = []
						if item != "NA" and item != "NaN" and item != "nan" and item != "Inf" and item != "inf":
							total_value[cluster][mycmp]["yes"].append(float(item))
						if not cluster in abundance:
							abundance[cluster] = {}
						abundance[cluster][mys] = item
						if not mycmp in meta_case_new:
							meta_case_new[mycmp] = {}
						meta_case_new[mycmp][mys] = ""
				if mycmp in meta_control:
					if mys in meta_control[mycmp]:
						if not cluster in total_value:
							total_value[cluster] = {}
						if not mycmp in total_value[cluster]:
							total_value[cluster][mycmp] = {}
						if not "no" in total_value[cluster][mycmp]:
							total_value[cluster][mycmp]["no"] = []
						if item != "NA" and item != "NaN" and item != "nan" and item != "Inf" and item != "inf":
							total_value[cluster][mycmp]["no"].append(float(item))
						if not cluster in abundance:
							abundance[cluster] = {}
						abundance[cluster][mys] = item
						if not mycmp in meta_control_new:
							meta_control_new[mycmp] = {}
						meta_control_new[mycmp][mys] = ""
			myindex = myindex + 1
	# foreah sample
	open_file.close()
	sys.stderr.write("Collect abundance info......done\n")

	# get fold change
	sys.stderr.write("Collect abundance fold change......starting\n")
	folds = {}
	for myclust in total_value.keys():
		for mycmp in total_value[myclust].keys():
			myyes = 0
			myno = 0
			logyes = 0
			logno = 0
			myyes_sd = 0
			myno_sd = 0
			if "yes" in total_value[myclust][mycmp]:
				if len(total_value[myclust][mycmp]["yes"]) > 0:
					#myyes = utilities.mean(total_value[myclust][mycmp]["yes"])
					myyes = 1.0 * sum(total_value[myclust][mycmp]["yes"]) / len(total_value[myclust][mycmp]["yes"])
					logyes = [math.log(k) for k in total_value[myclust][mycmp]["yes"]]
					logyes = 1.0 * sum(logyes) / len(total_value[myclust][mycmp]["yes"])
				if len(total_value[myclust][mycmp]["yes"]) > 1:
					myyes_sd = utilities.stddev(total_value[myclust][mycmp]["yes"], ddof=1)
			if "no" in total_value[myclust][mycmp]:
				if len(total_value[myclust][mycmp]["no"]) > 0:
					#myno = utilities.mean(total_value[myclust][mycmp]["no"])
					myno = 1.0 * sum(total_value[myclust][mycmp]["no"]) / len(total_value[myclust][mycmp]["no"])
					logno = [math.log(k) for k in total_value[myclust][mycmp]["no"]]
					logno = 1.0 * sum(logno) / len(total_value[myclust][mycmp]["no"])
				if len(total_value[myclust][mycmp]["no"]) > 1:
					myno_sd = utilities.stddev(total_value[myclust][mycmp]["no"], ddof=1)
			myyes_real = myyes
			myno_real = myno
			mean_case =  myyes
			mean_control = myno
			if myno == 0:
				if myyes == 0:
					myfold = "NaN"
				else:
					myfold = "Inf"
			else:
				myfold = math.log(myyes / (myno * 1.0))
			myvar = logyes - logno
			if myyes_sd == 0 and myno_sd == 0:
				# debug
				print("SD is zero!\t" + mycmp + "\t" + myclust + "\t" + str(myyes_real) + "\t" + str(myno_real))
				myeffect = "NaN"
			else:
				myeffect = (myyes_real - myno_real) / math.sqrt((myyes_sd * myyes_sd + myno_sd * myno_sd) / 2)
			
			# mean prevalence abundance
			#tmp1 = [x for x in total_value[myclust][mycmp]["yes"] if float(x) > float(config.abundance_detection_level)]
			tmp1 = [x for x in total_value[myclust][mycmp]["yes"] if float(x) != 0]
			tmp2 = [x for x in total_value[myclust][mycmp]["no"] if float(x) != 0]
			if len(tmp1) > 0:
				prevalent_mean_case = 1.0 * sum(tmp1)/len(tmp1)
			else:
				prevalent_mean_case = 0
			if len(tmp2) > 0:
				prevalent_mean_control = 1.0 * sum(tmp2)/len(tmp2)
			else:
				prevalent_mean_control = 0
			prevalent_tmp = []
			prevalent_tmp.extend(tmp1)
			prevalent_tmp.extend(tmp2)
			tmp = []
			tmp.extend(total_value[myclust][mycmp]["yes"])
			tmp.extend(total_value[myclust][mycmp]["no"])
			if len(tmp) == 0:
				print("Zero prevalence\t" + myclust)
				mymean = 0
				myvar = 0
			else:
				mymean = 1.0 * sum(tmp)/len(tmp)
			if len(prevalent_tmp) == 0:
				prevalent_mymean = 0
			else:
				prevalent_mymean = 1.0 * sum(prevalent_tmp)/len(prevalent_tmp)
			if not myclust in folds:
				folds[myclust] = {}
			folds[myclust][mycmp] = str(myfold) + "\t" + str(myeffect) + "\t" + str(myvar) + "\t" + str(mymean) + "\t" + str(mean_case) + "\t" + str(mean_control) + "\t" + str(prevalent_mymean) + "\t" + str(prevalent_mean_case) + "\t" + str(prevalent_mean_control)
	# foreach cluster
	total_value = {}
	tmp1 = []
	tmp2 = []
	tmp = []
	sample = {}
	sys.stderr.write("Collect abundance fold change......done\n")
	return folds, abundance, meta_case_new, meta_control_new

# function collect_abundance_info_fold


# ==============================================================
# collect prevalence
# ==============================================================
def collect_prevalence_info (abundance_file, DA_cluster, meta_case, meta_control):  # summary_peptide_family_abundance.tsv
	pre_meta_case = {}
	pre_meta_control = {}
	samples = {}
	open_file = open(abundance_file, "r")
	sys.stderr.write("Collect prevalence info......starting\n")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 1
	while myindex < len(info):
		item = info[myindex]
		samples[myindex] = item
		myindex = myindex + 1
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		cluster = info[0]
		if not cluster in DA_cluster:
			continue
		cmp_types = DA_cluster[cluster].keys()
		myindex = 1
		while myindex < len(info):
			item = info[myindex]
			mys = samples[myindex]
			for mycmp in cmp_types:
				if mycmp in meta_case:
					if mys in meta_case[mycmp]:
						if item != "NA" and item != "NaN" and item != "nan" and item != "Inf" and item != "inf":
							count = float(item)
							#if count > float(config.abundance_detection_level):
							if count != 0:
								if not cluster in pre_meta_case:
									pre_meta_case[cluster] = {}
								if not mycmp in pre_meta_case[cluster]:
									pre_meta_case[cluster][mycmp] = {}
								pre_meta_case[cluster][mycmp][mys] = ""
				if mycmp in meta_control:
					if mys in meta_control[mycmp]:
						if item != "NA" and item != "NaN" and item != "nan" and item != "Inf" and item != "inf":
							count = float(item)
							#if count > float(config.abundance_detection_level):
							if count != 0:
								if not cluster in pre_meta_control:
									pre_meta_control[cluster] = {}
								if not mycmp in pre_meta_control[cluster]:
									pre_meta_control[cluster][mycmp] = {}
								pre_meta_control[cluster][mycmp][mys] = ""
			myindex = myindex + 1
	# foreach family
	open_file.close()
	sys.stderr.write("Collect prevalence info......done\n")

	# get prevalence
	prevalence = {}
	sys.stderr.write("Calculate prevalence info......starting\n")
	for mycmp in sorted(meta_control.keys()):
		meta_yes_num = len(meta_case[mycmp].keys())
		meta_no_num = len(meta_control[mycmp].keys())
		total_num = meta_yes_num + meta_no_num
		#print(">>" + mycmp + "\t" + str(meta_yes_num) + "\t" + str(meta_no_num) + "\t" + str(total_num))

		for myclust in DA_cluster.keys():
			mymeta_yes = 0
			mymeta_no = 0
			if myclust in pre_meta_case:
				if mycmp in pre_meta_case[myclust]:
					mymeta_yes = len(pre_meta_case[myclust][mycmp].keys())
			if myclust in pre_meta_control:
				if mycmp in pre_meta_control[myclust]:
					mymeta_no = len(pre_meta_control[myclust][mycmp].keys())
			#print(myclust + "\t" + str(mymeta_yes) + "\t" + str(mymeta_no))
			pre_yes = 1.0 * mymeta_yes / meta_yes_num
			pre_no = 1.0 * mymeta_no / meta_no_num
			pre = 1.0 * (mymeta_yes + mymeta_no) / total_num
			if not myclust in prevalence:
				prevalence[myclust] = {}
			prevalence[myclust][mycmp] = str(pre) + "\t" + str(pre_yes) + "\t" + str(pre_no)
	# foreach cluster
	sys.stderr.write("Calculate prevalence info......done\n")
	pre_meta_case = {}
	pre_meta_control = {}
	samples = {}
	return prevalence

# function get_prevalence



#==============================================================
# Output info 
#==============================================================
def output_DA_info (stats, outfile):
	outfile = re.sub(".tsv", ".stat.tsv", outfile)
	open_file = open(outfile, "w")
	open_file.write(utilities.PROTEIN_FAMILY_ID + "\tcmp_type\tprevalence\tcoef\tstderr\tpval\tqval\n")
	for myclust in sorted(stats.keys()):
		for mycmp in sorted(stats[myclust].keys()):
			open_file.write(myclust + "\t" + mycmp + "\t" + stats[myclust][mycmp] + "\n")
	# foreach cluster
	open_file.close()
# output_info

# output abundance
def output_abundance_info (folds, abundance, meta_control, meta_case, outfile):
	'''
	for mycmp in sorted(meta_control.keys()):
		title = "ID"
		title = title + "\t" + "\t".join(meta_case[mycmp].keys())
		title = title + "\t" + "\t".join(meta_control[mycmp].keys())
		outfile1 = re.sub(".tsv", "." + mycmp + ".abundance.tsv", outfile)
		open_file = open(outfile1, "w")
		open_file.write(title + "\n")
		for myclust in sorted(abundance.keys()):
			mystr = myclust
			for mys in sorted(meta_case[mycmp].keys()):
				if mys in abundance[myclust]:
					mystr = mystr + "\t" + str(abundance[myclust][mys])
				else:
					# debug
					#print("Contrast sample abundance not exist!\t" + myclust + "\t" + mys)
					mystr = mystr + "\t0"
			for mys in sorted(meta_control[mycmp].keys()):
				if mys in abundance[myclust]:
					mystr = mystr + "\t" + str(abundance[myclust][mys])
				else:
					# debug
					#print("Ref sample abundance not exist!\t" + myclust + "\t" + mys)
					mystr = mystr + "\t0"
			open_file.write(mystr + "\n")
		open_file.close()
	# foreach cmp type
	'''

	data_file = re.sub(".tsv", ".fold.tsv", outfile)
	title = utilities.PROTEIN_FAMILY_ID + "\tcmp_type\tlog(FC)\tCohen's_d\tmean(log)\tmean_abundance\tmean_abundance_case\tmean_abundance_control\tmean_prevalent_abundance\tmean_prevalent_abundance_case\tmean_prevalent_abundance_control"
	open_file = open(data_file, "w")
	open_file.write(title + "\n")
	for myid in sorted(folds.keys()):
		for mycmp in sorted(folds[myid].keys()):
			open_file.write(myid + "\t" + mycmp + "\t" + str(folds[myid][mycmp]) + "\n")
	open_file.close()
# output_abundance_info

# output prevalence
def output_prevalence_info (prevalence, outfile):
	outfile = re.sub(".tsv", ".prevalence.tsv", outfile)
	title = utilities.PROTEIN_FAMILY_ID + "\tcmp_type\tprevalence\tprevalence_case\tprevalence_control"
	open_file = open(outfile, "w")
	open_file.write(title + "\n")
	for myclust in sorted(prevalence.keys()):
		for mycmp in sorted(prevalence[myclust].keys()):
			open_file.write(myclust + "\t" + mycmp + "\t" + prevalence[myclust][mycmp] + "\n")
	# foreach cluster
	open_file.close()
# output_prevalence_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start maaslin2_collection.py -i " + values.maaslin2_results+ " ####\n")
	
	
	### collect abundance info ###
	sys.stderr.write("Get stat info......starting\n")
	DA_cluster = collect_DA_info (values.maaslin2_results)
	if values.metadata == "none":
		meta_case, meta_control = collect_meta_info (config.metadata)
	else:
		meta_case, meta_control = collect_meta_info (values.metadata)
	folds, abundance, meta_case, meta_control = collect_abundance_info (values.smoothed_abundance, DA_cluster, meta_case, meta_control)
	prevalence = collect_prevalence_info (values.abundance_table, DA_cluster, meta_case, meta_control)
	sys.stderr.write("Get stat info......done\n")

	### output info ###
	sys.stderr.write("Output stat file......starting\n")
	output_DA_info (DA_cluster, values.output)
	output_abundance_info(folds, abundance, meta_control, meta_case, values.output)
	output_prevalence_info(prevalence, values.output)
	sys.stderr.write("Output stat file......done\n")
	

	sys.stderr.write("### Finish maaslin2_collection.py ####\n\n\n")

# end: main

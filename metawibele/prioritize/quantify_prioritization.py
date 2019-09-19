#!/usr/bin/env python

"""
MeteWIBELE: quantify_prioritization module
1) Define quantitative criteria to prioritize the importance of protein families
2) Prioritize the importance of protein families using unsupervised or supervised approaches

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
import argparse
import subprocess
import tempfile
import re
import logging
import numpy
import scipy.stats
import pandas as pd
from collections import namedtuple
from operator import attrgetter, itemgetter

# Try to load one of the MetaWIBELE modules to check the installation
try:
	import config
	import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# name global logging instance
logger = logging.getLogger(__name__)


def parse_arguments():
	"""
	Parse the arguments from the user
	"""
	parser = argparse.ArgumentParser(
		description = "MetaWIBELE-prioritize: prioritize importance of protein families based on annotation and phenotype properties\n",
		formatter_class = argparse.RawTextHelpFormatter,
		prog = "prioritization")
	parser.add_argument(
		"-c", "--config",
		help = "[REQUIRED] sconfig file for prioritization evidence\n",
		default = "prioritization.cfg",
		required=True)
	parser.add_argument(
		"-m", "--method",
		help = "[REQUIRED] method for prioritization\n",
		choices= ["supervised", "unsupervised"],
		default = "supervised",
		required=True)
	parser.add_argument(
		"-w", "--weight",
		help = "[REQUIRED] method for weighting: "
		       "[uniform] specify uniform weight for each evidence; "
		       "[correlated] specify weigh based on the pairwise correlation between evidence items;"
		       "[fixed] specify weigh manually in the config file\n",
		choices= ["uniform", "correlated", "fixed"],
		default = "uniform",
		required=True)
	parser.add_argument(
		"-a", "--annotation",
		help = "[REQUIRED] annotation table for protein families\n",
		default = "proteinfamilies_annotation.tsv",
		required=True)
	parser.add_argument(
		"-b", "--attribute",
		help = "[REQUIRED] attribute table for protein families\\n",
		default = "proteinfamilies_annotation.attribute.tsv",
		required=True)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] writing directory for output files\n",
		default = "prioritization",
		required=True)

	return parser.parse_args()


def read_config_file (conf_file, method):
	"""
	Collect config info for prioritization
	Input: config filename
	Output: evidence_conf = {DNA_prevalence:1, DNA_abundance:1, ...}
	"""

	logger.debug('read_config_file')

	config_items = config.read_user_edit_config_file(conf_file)
	ann_conf = {}
	attr_conf = {}

	if method == "unsupervised":
		if "unsupervised" in config_items:
			for name in config_items["unsupervised"].keys():
				myvalue = config_items["unsupervised"][name]
				try:
					float(myvalue)
					if float(myvalue) == 0:
						continue
				except ValueError:
					logger.debug('Not numberic values for the config item')
				if re.search("__", name):
					name = re.sub("-", "_", name)
					name = re.sub("\.", "_", name)
					attr_conf[name] = myvalue
				else:
					name = re.sub("-", "_", name)
					name = re.sub("\.", "_", name)
					ann_conf[name] = myvalue

	if method == "supervised":
		if "supervised" in config_items:
			for name in config_items["supervised"].keys():
				myvalue = config_items["supervised"][name]
				try:
					float(myvalue)
					if float(myvalue) == 0:
						continue
				except ValueError:
					logger.debug('Not numberic values for the config item')
				if re.search("__", name):
					name = re.sub("-", "_", name)
					name = re.sub("\.", "_", name)
					attr_conf[name] = myvalue
				else:
					name = re.sub("-", "_", name)
					name = re.sub("\.", "_", name)
					ann_conf[name] = myvalue
	
	return ann_conf, attr_conf


def read_attribute_file (attr_file, attr_conf):
	"""
	Collect annotation evidence for protein families used for prioritization
	Input: filename of the characterization file
	Output: ann = {Cluster_XYZ: {qvalue:0.001, coef:-0.3, ...}, ...}
	"""
	
	required = {}
	annotation = {}
	split = {}
	flags = {}
	titles = {}
	open_file = open(attr_file, "r")
	line = open_file.readline()
	line = re.sub("\n$", "", line)
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = re.sub("\n$", "", line)
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[titles["AID"]]
		myclust, mytype = myid.split("__")
		myid = myclust
		mykey = info[titles["key"]]
		mytype_new = mytype + "__" + mykey
		mytype_new = re.sub("-", "_", mytype_new)
		mytype_new = re.sub("\.", "_", mytype_new)
		myvalue = info[titles["value"]]
		if mykey == "cmp_type":
			flags[myid] = myvalue
		if not mytype_new.lower() in attr_conf:
			continue
		if attr_conf[mytype_new.lower()] == "required":
			required[mytype_new] = ""
		if re.search("MaAsLin2", mytype) and myid in flags:
			myclust = myid + "|" + flags[myid]
			if not myid in split:
				split[myid] = {}
			split[myid][myclust] = ""
		if myvalue == "NA" or myvalue == "NaN" or myvalue == "nan" or myvalue == "Nan":
			continue
		if not myclust in annotation:
			annotation[myclust] = {}
		annotation[myclust][mytype_new] = myvalue
	# foreach line
	open_file.close()

	return annotation, split, required


def read_annotation_file (ann_file, ann_conf):
	"""
	Collect annotation evidence for protein families used for prioritization
	Input: filename of the characterization file
	Output: ann = {Cluster_XYZ: {prevalence:0.001, abundance:0.3, ...}, ...}
	"""
	logger.debug('read_annotation_file')

	required = {}
	annotation = {}
	titles = {}
	open_file = open(ann_file, "r")
	line = open_file.readline()
	line = re.sub("\n$", "", line)
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = re.sub("\n$", "", line)
		if not len(line):
			continue
		info = line.split("\t")
		myclust = info[titles[utilities.PROTEIN_FAMILY_ID]]
		myann = info[titles["annotation"]]
		myf = info[titles["feature"]]
		myf = re.sub("-", "_", myf)
		myf = re.sub("\.", "_", myf)
		if myann == "NA" or myann == "NaN" or myann == "nan" or myann == "Nan":
			continue
		if myf.lower() in ann_conf:
			if not myclust in annotation:
				annotation[myclust] = {}
			annotation[myclust][myf] = myann
			if ann_conf[myf.lower()] == "required":
				required[myf] = ""
	# foreach line
	open_file.close()

	return annotation, required


def combine_annotation (annotation, split, required, total_ann, ann_types, required_types):
	"""
	Combine annotation information of protein families for prioritization
	Input: ann = {Cluster_XYZ: {prevalence:0.001, abundance:0.3, ...}, ...}
			attr = {Cluster_XYZ: {prevalence:0.001, abundance:0.3, ...}, ...}
			split = {Cluster_XYZ:{Cluster_XYZ|A, Cluster_XYZ|B, ...}, ...}
	Output: total = {Cluster_XYZ: {prevalence:0.001, abundance:0.3, ...}, ...}
	"""
	logger.debug('combine_annotation')

	for myid in annotation.keys():
		if myid in split:
			for myid_new in split[myid].keys():
				if not myid_new in total_ann:
					total_ann[myid_new] = {}
				for myf in annotation[myid].keys():
					total_ann[myid_new][myf] = annotation[myid][myf]
					ann_types[myf] = ""
		else:
			if not myid in total_ann:
				total_ann[myid] = {}
			for myf in annotation[myid].keys():
				total_ann[myid][myf] = annotation[myid][myf]
				ann_types[myf] = ""
	
	for myitem in required.keys():
		required_types[myitem] = ""


def check_annotation (annotation, required_types):
	"""
	Select clusters with required annotation types
	Input: ann = {Cluster_XYZ: {prevalence:0.001, abundance:0.3, ...}, ...}
	Output: ann_new = {Cluster_abc: {prevalence:0.001, abundance:0.3, ...}, ...}
	"""

	# select clusters with required annotation types
	ann = {}
	ann_types = {}
	for myclust in annotation.keys():
		myflag = 0
		for myitem in required_types.keys():
			if not myitem in annotation[myclust]:
				#print("No required type\t" + myitem + "\t" + myclust)
				myflag = 1
				break
		if myflag == 0:
			if not myclust in ann:
				ann[myclust] = {}
			for myitem in annotation[myclust].keys():
				 ann[myclust][myitem] = annotation[myclust][myitem]
				 ann_types[myitem] = ""

	return ann, ann_types


def combine_evidence (ann, ann_types):
	"""
	Combine prioritization evidence for protein families
	Input: ann = {Cluster_XYZ: {'qvalue':0.001, 'coef':-0.3, ...}, ...}
		   ann_types = {'qvalue', 'coef', ...}
	Output: evidence_dm = {Cluster_XYZ: {'qvalue':0.001, 'coef':-0.3, 'annotation':3, ...}, ...}
	"""

	logger.debug('combine_evidence')

	evidence_row = sorted(ann_types.keys())
	metawibele_row = []
	for item in evidence_row:
		metawibele_row.append(item + "__value")
		metawibele_row.append(item + "__percentile")
	evidence_table_row = namedtuple("evidence_table_row", evidence_row, verbose=False, rename=False)
	evidence_table = pd.DataFrame(index=sorted(ann.keys()), columns=evidence_table_row._fields)

	# build data frame
	for item in evidence_row:
		myvalue = []
		for myclust in sorted(ann.keys()):
			if item in ann[myclust]:
				myvalue.append(ann[myclust][item])
			else:
				# debug
				#print("No item!\t" + myclust + "\t" + item)
				myvalue.append("NaN")
		# foreach cluster
		evidence_table[item] = myvalue
	# foreach evidence

	return evidence_table, evidence_row, metawibele_row


def get_correlated_weight (evidence_table):
		"""
		Calculate the pairwise correlation between evidence items and return weight table
		Input:  evidence_table = {family: {'abundance': abundance, 'prevalence': prevalence}}
		Output: weight_conf = {'abundance': 0.5, 'prevalence': 0.5, ...}
		"""

		df = evidence_table
		df = df.apply(pd.to_numeric, errors='coerce')
		weight_conf = {}
		df_corr = df.corr(method="spearman")
		df_corr = abs(df_corr)
		df_corr['weight'] = 1.0 / df_corr.sum(skipna=True)
		for index, row in df_corr.iterrows():
			weight_conf[index] = row.weight
			print(index + "\t" + str(row.weight))
		
		return weight_conf


def get_uniform_weight (ann_types):
	"""
	Calculate the uniform weight and return weight table
	Input:  evidence_table = {family: {'abundance': abundance, 'prevalence': prevalence}r
	Output: weight_conf = {'abundance': 0.5, 'prevalence': 0.5, ...}
	"""

	weight_conf = {}
	myweight = 1.0 / len(ann_types.keys())
	for mytype in ann_types.keys():
		weight_conf[mytype] = myweight
		print(mytype + "\t" + str(myweight))
	
	return weight_conf


def get_fixed_weight (ann_types, ann_conf, attr_conf):
	"""
	Calculate the uniform weight and return weight table
	Input:  evidence_table = {family: {'abundance': abundance, 'prevalence': prevalence}}
	Output: weight_conf = {'abundance': 0.5, 'prevalence': 0.5, ...}
	"""
	
	weight_conf = {}
	for mytype in ann_types.keys():
		if mytype.lower() in ann_conf:
			weight_conf[mytype] = ann_conf[mytype.lower()]
			# debug
			print(mytype + "\t" + str(ann_conf[mytype.lower()]))
		if mytype.lower() in attr_conf:
			weight_conf[mytype] = attr_conf[mytype.lower()]
			print(mytype + "\t" + str(attr_conf[mytype.lower()]))
	
	return weight_conf


def weighted_harmonic_mean (summary_table, evidence, weight_conf, score_name):
	"""
	Calculate the weighted harmonic mean
	Input:  summary_table = {family: {'abundance': 0.5, 'prevalence': 0.8}, ...}
			evidence = ['abundance', 'prevalence', ...]
			weight_conf = {'abundance': 0.5, 'prevalence': 0.5, ...}
	Output: summary_table = {family: {'score_name': 0.9, 'abundance_value': 0.5, 'abundance_percentile':0.9,...},...}
	"""

	# Weighted Harmonic mean
	total_weight = 0
	mytype = evidence[0]
	mykey = mytype + "__percentile"
	myw = float(weight_conf[mytype])
	total_weight = total_weight + myw
	myscore = myw / summary_table[mykey]
	for mytype in evidence[1:]:
		mykey = mytype + "__percentile"
		if mytype in weight_conf:
			myw = float(weight_conf[mytype])
			total_weight = total_weight + myw
			myscore = myscore + myw / summary_table[mykey]
	summary_table[score_name] = float(total_weight) / myscore
	

def get_rank_score (evidence_table, evidence_row, metawibele_row, weight_conf):
		"""
		Return the data frame of protein families with their annotation, percentiles, and MetaWIBELE score
		Input:  evidence_table = {family: {'abundance': 0.5, 'prevalence': 0.8}}
				beta = parameter value
		Output: summary_table = {family: {'abundance_value': 0.5, 'abundance_percentiles': 0.9,...},...}
		"""

		logger.debug('get_rank_score')

		# create a data frame
		metawibele_table_row = namedtuple("metawibele_table_row", metawibele_row, verbose=False, rename=False)
		summary_table = pd.DataFrame(index=evidence_table.index, columns=metawibele_table_row._fields)

		# calculate percentile
		rank_name = []
		for mytype in evidence_row:
			summary_table[mytype + "__value"] = evidence_table[mytype]
			summary_table[mytype + "__percentile"] = scipy.stats.rankdata(pd.to_numeric(summary_table[mytype + "__value"], errors='coerce'), method='average')
			if re.search(config.effect_size, mytype):
				summary_table[mytype + "__percentile"] = scipy.stats.rankdata(abs(pd.to_numeric(summary_table[mytype + "__value"], errors='coerce')), method='average')
			else:
				if re.search("qvalue", mytype) or re.search("q-value", mytype) or re.search("pvalue", mytype) or re.search("p-value", mytype):
					summary_table[mytype + "__percentile"] = scipy.stats.rankdata(-pd.to_numeric(summary_table[mytype + "__value"], errors='coerce'), method='average')
			summary_table[mytype + "__percentile"] = summary_table[mytype + "__percentile"] / summary_table[mytype + "__percentile"].max()
			rank_name.append(mytype + "__percentile")

		# calculate MetaWIBELE score
		weighted_harmonic_mean (summary_table, evidence_row, weight_conf, "ranking_score")
		#summary_table["ranking_score"] = summary_table[rank_name].min(axis=1)
		summary_rank = summary_table[rank_name]

		return summary_table, summary_rank


def prioritize_families (summary_table, score_column, ann_conf):
	"""
	Return important protein families based on MetaWIBELE score
	Input: summary_table = {family: {'abundance': mean abundance, 'prevalence': prevalence}}
			beta = parameter value
	Output: imp_families = {family: {'abundance': mean abundance, 'prevalence': prevalence}}
	"""
	logger.debug('prioritize_families')

	pri_beta = ann_conf["tshld_priority"]
	pri_score = ann_conf["tshld_priority_score"]
	#metawibele_score = config.tshld_score  # 1/((1/(beta*tshld_prev)) + (1/((1-beta)*tshld_abund)))

	# get important family based on their MetaWIBELE score
	try:
		# for  pandas >= 0.17.0
		#imp_families = imp_families.sort_values(by=score_column, ascending=False)
		summary_table = summary_table.sort_values(by=score_column, ascending=False)
	except:
		#imp_families = imp_families.sort(score_column, ascending=False)
		summary_table = summary_table.sort_values(score_column, ascending=False)
	if pri_score != "none":
		imp_families = summary_table[summary_table[score_column] > float(pri_score)]
		# debug
		print("Specified threshold of priority score: " + str(pri_score))
	else:
		mytop_num = int(summary_table.shape[0] * float(pri_beta))
		imp_families = summary_table.head(mytop_num)
		print("Specified threshold of priority: " + str(pri_beta))
	
	return summary_table, imp_families


def write_results (summary_table, split, out_file):
	"""
	Write the prevalence, abundance, rank score information for protein families in text file
	Input: summary_table = {family: {'mean_abundance': mean abundance, ...}}
			out_file = output_filename
	Output: Writes the family dictionary to the output_filename
	"""

	logger.debug('write_prioritization_results')

	keys = summary_table.columns.values.tolist()
	foo = open(out_file, 'w')
	foo.write(utilities.PROTEIN_FAMILY_ID + "\t" + "\t".join(keys) + "\n")
	#if len(split.keys()) > 0:
	#	foo.write(utilities.PROTEIN_FAMILY_ID + "\t" + "\t".join(keys) + "\n")
	#else:	
	#	foo.write(utilities.PROTEIN_FAMILY_ID + "\t" + "\t".join(keys) + "\n")
	foo.close()
	summary_table.to_csv(out_file, mode='a', sep='\t', header=False)


def main():
	args_value = parse_arguments()
	myout = config.basename + "_" + args_value.method + "_prioritization"

	### calculate Ranking score ###
	if config.verbose == 'DEBUG':
		print ("--- Collecting annotations for protein families ---")
	ann_conf, attr_conf = read_config_file (args_value.config, args_value.method)
	attribute, split, required_attr = read_attribute_file (args_value.attribute, attr_conf)
	annotation, required_ann = read_annotation_file (args_value.annotation, ann_conf)
	ann = {}
	ann_types = {}
	required_types = {}
	combine_annotation (attribute, split, required_attr, ann, ann_types, required_types)
	combine_annotation (annotation, split, required_ann, ann, ann_types, required_types)
	ann_new, ann_types_new = check_annotation (ann, required_types)
	evidence_table, evidence_row, metawibele_row = combine_evidence (ann_new, ann_types_new)

	if config.verbose == 'DEBUG':
		print ("--- MetaWIBELE evidence table are written to the output ---")
	if not os.path.exists(args_value.output):
		os.system("mkdir -p " + args_value.output)
	metawibele_output_file = args_value.output + '/' + myout + '.evidence.tsv'
	write_results (evidence_table, split, metawibele_output_file)

	if config.verbose == 'DEBUG':
		print ("--- Calculate Ranking score for protein families ---")
	if args_value.weight == "fixed":
		weight_conf = get_fixed_weight (ann_types, ann_conf, attr_conf)
	if args_value.weight == "uniform":
		weight_conf = get_uniform_weight (ann_types)
	if args_value.weight == "correlated":
		weight_conf = get_correlated_weight (evidence_table)
	summary_table, rank_table = get_rank_score (evidence_table, evidence_row, metawibele_row, weight_conf)

	### get important families ###
	if config.verbose == 'DEBUG':
		print ("--- Get prioritized families based on MetaWIBELE score ---")
	summary_table, imp_family = prioritize_families (summary_table, "ranking_score", ann_conf)
	metawibele_output_file = args_value.output + '/' + myout + '.rank.tsv'
	write_results(summary_table, split, metawibele_output_file)
	metawibele_output_file = args_value.output + '/' + myout + '.priority.tsv'
	write_results (imp_family, split, metawibele_output_file)

	if config.verbose == 'DEBUG':
		print ("--- The prioritization output is written in %s ..." % (args_value.output))
		print ("--- Prioritization process is successfully completed ---")


if __name__ == '__main__':
	main()

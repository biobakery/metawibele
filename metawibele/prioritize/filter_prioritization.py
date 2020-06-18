#!/usr/bin/env python

"""
MeteWIBELE: filter_prioritization.py module
1) Filter subsets using binary annotation (sequence-based annotation, interested genes, etc) as prioritization factor 
2) Narrow down prioritized candidates

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
	from metawibele import config
	from metawibele import utilities
	from metawibele.common import utils
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
		description = "MetaWIBELE-prioritize: prioritize importance of protein families based on binary annotation\n",
		formatter_class = argparse.RawTextHelpFormatter,
		prog = "filter_prioritization.py")
	parser.add_argument(
		"-c", "--config",
		help = "[REQUIRED]config file for prioritization filtering\n",
		default = "prioritization.cfg",
		required=True)
	parser.add_argument(
		"-f", "--function",
		help = "[OPTIONAL] specify file for interested functions filtering\n",
		default = "none")
	parser.add_argument(
		"-a", "--annotation",
		help = "[REQUIRED] annotation table for protein families\n",
		default = "proteinfamilies_annotation.tsv",
		required=True)
	parser.add_argument(
		"-p", "--priority",
		help = "[REQUIRED] prioritization file\n",
		default = "metawibele_supervised_prioritization.priority.tsv",
		required=True)
	parser.add_argument(
		"-o", "--outfile",
		help = "[REQUIRED] the name of output file\n",
		default = "metawibele_supervised_prioritization.priority-select.tsv",
		required=True)

	return parser.parse_args()


def read_config_file (conf_file):
	"""
	Collect config info for prioritization filtering
	Input: config filename
	Output: evidence_conf = {signaling:1, interaction:1, extracellular:1 ...}
	"""

	print('read_config_file')

	config_items = config.read_user_edit_config_file(conf_file)
	required_conf = {}
	optional_conf = {}
	specific_annotation = {}
	specific_cluster = {}
	values = ["required", "optional", "none"]

	if "filtering" in config_items:
		for name in config_items["filtering"].keys():
			myvalue = config_items["filtering"][name]
			if name == "vignettes" or name == "clusters":
				if myvalue.lower() == "none":
					continue
				print("Required filtering item: " + name + "\t" + myvalue)
				required_conf[name] = myvalue
				if name == "vignettes":
					specific_annotation[myvalue] = ""
				if name == "clusters":
					specific_cluster[myvalue] = ""
			else:
				if not myvalue in values:
					print("Please use valid value for the config item " + name + ": e.g. required | optional | none")
					continue
				if myvalue.lower() == "none":
					continue
				if myvalue.lower()  == "required":
					print("Required filtering item: " + name + "\t" + myvalue)
					required_conf[name] = myvalue
				if myvalue.lower() == "optional": 
					print("Optional filtering item: " + name + "\t" + myvalue)
					optional_conf[name] = myvalue

	return required_conf, optional_conf, specific_annotation, specific_cluster


def read_vignettes_file (vignettes_file, specific_annotation):
	"""
	Collect proteins with specific annotations for prioritization
	Input: vignettes filename
	Output: vignettes = [ann1, ann2, ann3, ..]
	"""
	print('read_vignettes_file')

	vignettes = {}
	titles = {}
	for line in utils.gzip_bzip2_biom_open_readlines (vignettes_file):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^#", line):
			continue
		info = line.split("\t")
		if re.search("^type", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		mytype = info[titles["type"]]
		if not mytype.lower() in specific_annotation and not mytype in specific_annotation:
			continue
		if not "annotation" in titles:
			# debug
			print("No annotation info!\t" + line)
			continue
		myid = info[titles["annotation"]]
		vignettes[myid] = ""
	# foreach line

	return vignettes


def read_cluster_file (specific_cluster):
	"""
	Collect proteins with specific clusters for prioritization
	Input: cluster filename
	Output: clusters = [ann1, ann2, ann3, ..]
	"""

	print('read_cluster_file')

	cluster = {}
	titles = {}
	for myfile in specific_cluster.keys():
		if not myfile.exists():
			# debug
			print("File not exist!\t" + myfile)
			continue
		open_file = open(myfile, "r")
		line = open_file.readline()
		line = re.sub("\n$", "", line)
		info = line.split("\t")
		for item in info:
			titles[item] = info.index(item)
		for line in open_file:
			line = re.sub("\n$", "", line)
			if not len(line):
				continue
			if re.search("^#", line):
				continue
			info = line.split("\t")
			myid = info[titles[utilities.PROTEIN_FAMILY_ID]]
			cluster[myid] = ""
		# foreach line
		open_file.close()
	# foreach file

	return cluster



def read_annotation_file (ann_file, required_conf, optional_conf, specific_ann, specific_cluster):
	"""
	Collect annotation evidence for protein families
	Input: filename of the characterization file
	Output: ann_score = {Cluster_XYZ: 3, ...}
	"""

	print('read_annotation_file')

	ann = {}
	ann_types = {}
	cluster = {}
	titles = {}
	open_file = open(ann_file, "r")
	for line in open_file:
		line = re.sub("\n$", "", line)
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^familyID", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myid = info[titles["familyID"]]
		myann = info[titles["annotation"]]
		myf = info[titles["feature"]]
		if myann == "NA":
			continue
		try:
			float(myann)
			if float(myann) == 0:
				continue
		except ValueError:
			#print('Not numberic values for annotation item' + "\t" + myann)
			pass
		if myf.lower() in required_conf or myf.lower() in optional_conf:
			if not myid in ann:
				ann[myid] = {}
			if not myf in ann[myid]:
				ann[myid][myf] = myann
		if len(specific_ann.keys()) > 0:
			tmp = myann.split(";")
			for item in tmp:
				if item in specific_ann:
					if not myid in ann:
						ann[myid] = {}
					if not "vignettes" in ann[myid]:
						ann[myid]["vignettes"] = 1
		if len(specific_cluster.keys()) > 0:
			if myid in specific_cluster:
				if not myid in ann:
					ann[myid] = {}
				if not "clusters" in ann[myid]:
					ann[myid]["clusters"] = 1
	# foreach line
	open_file.close()
	
	# check required and optional annotation
	cluster1 = {}
	for myid in ann.keys():
		mynum = 0
		for mytype in ann[myid].keys():
			if mytype.lower() in required_conf.keys():
				mynum = mynum + 1
			ann_types[mytype] = ""
		if len(required_conf.keys()) == mynum:
			cluster1[myid] = ""

	cluster2 = {}
	for myid in ann.keys():
		flag = 0
		for mytype in ann[myid].keys():
			if mytype.lower() in optional_conf.keys():
				flag = 1
			ann_types[mytype] = ""
		if len(optional_conf.keys()) == 0:
			flag = 1
		if flag == 1:
			cluster2[myid] = ""

	for myid in cluster1.keys():
		if myid in cluster2:
			cluster[myid] = ""

	return cluster


def filter_prioritization (priority_file, cluster, outfile):
	"""
	Filter prioritized protein families based on guiding evidence
	Input: prioritized file
			specific_families = {cluster_A, cluster_B, ...}
	Output: prioritized list after filtering
	"""

	print('filter_prioritization')

	open_file = open(priority_file, "r")
	open_out = open(outfile, "w")
	line = open_file.readline()
	open_out.write(line)
	titles = {}
	line = re.sub("\n$", "", line)
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = re.sub("\n$", "", line)
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[titles[utilities.PROTEIN_FAMILY_ID]]
		tmp = myid.split("|")
		myid = tmp[0]
		if myid in cluster:
			open_out.write(line + "\n")
	# foreach line
	open_file.close()
	open_out.close()


def main():
	args_value = parse_arguments()
	if args_value.function == "none":
		vignettes_database = config.vignettes_database
	else:
		vignettes_database = args_value.function

	if config.verbose == 'DEBUG':
		print ("--- Collect config information ---")
	required_conf, optional_conf, specific_annotation, specific_cluster = read_config_file (args_value.config)
	spe_ann = read_vignettes_file (vignettes_database, specific_annotation)
	spe_cluster = read_cluster_file (specific_cluster)
	clusters = read_annotation_file (args_value.annotation, required_conf, optional_conf, spe_ann, spe_cluster)

	if config.verbose == 'DEBUG':
		print ("--- Filter for subset of prioritized protein families ---")
	filter_prioritization (args_value.priority, clusters, args_value.outfile)


	if config.verbose == 'DEBUG':
		print ("--- The filter_prioritization output is written in %s ..." % (args_value.outfile))
		print ("--- Prioritization-filter process is successfully completed ---")


if __name__ == '__main__':
	main()

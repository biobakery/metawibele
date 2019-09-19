#!/usr/bin/env python
##########################################################################
# Function: Combine abundance details to specific features 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 08/29/2019
##########################################################################
import sys
import os
import os.path
import re
import argparse
import math

from metawibele import utilities

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', help='input DNA abundance file', required=True)
	parser.add_argument('-b', help='input RNA abundance file', required=True)
	parser.add_argument('-c', help='input RNA ratio file', required=True)
	parser.add_argument('-i', help='input features file', required=True)
	parser.add_argument('-o', help='output file', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect abundance 
#==============================================================
def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data)/n # in Python 2 use sum(data)/float(n)

def extract_abundance_info (abu_file):
	abundance = {}
	titles = {}
	open_file = open(abu_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_FAMILY_ID, line) or re.search("^" + utilities.PROTEIN_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myid = info[titles["familyID"]] + "\t" + info[titles["cmp_type"]]
		mypre = info[titles["prevalence"]]
		mypre_case = info[titles["prevalence_case"]]
		mypre_con = info[titles["prevalence_control"]]
		myfold = info[titles["foldChange"]]
		myabu_case = info[titles["mean_abundance_case"]]
		myabu_con = info[titles["mean_abundance_control"]]
		abundance[myid] = mypre_case + "\t" + mypre_con + "\t" + myabu_case + "\t" + myabu_con + "\t" + myfold
	# foreach line
	open_file.close()

	return abundance
	

#==============================================================
# combine information 
#==============================================================
def combine_info (dna_abu, rna_abu, ratio_abu, feature_file, outfile):
	titles = {}
	open_file = open(feature_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	title = line 
	title = re.sub("foldChange\tvalueChange\teffectSize\tabundance\t", "", title)
	title = re.sub("DNA_abundance\tDNA_prevalence", "DNA_prevalence_case\tDNA_prevalence_control\tDNA_mean_abundance_case\tDNA_mean_abundance_control\tDNA_foldChange",title)
	title = re.sub("RNA_abundance\tRNA_prevalence", "RNA_prevalence_case\tRNA_prevalence_control\tRNA_mean_abundance_case\tRNA_mean_abundance_control\tRNA_foldChange",title)
	title = re.sub("RNA-ratio_abundance\tRNA-ratio_prevalence", "RNA-ratio_prevalence_case\tRNA-ratio_prevalence_control\tRNA-ratio_mean_abundance_case\tRNA-ratio_mean_abundance_control\tRNA-ratio_foldChange",title)
	open_out = open(outfile, "w")
	open_out.write(title + "\n")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[titles["familyID"]] + "\t" + info[titles["cmp_type"]]
		mydna = "NA\tNA\tNA\tNA\tNA"
		myrna = "NA\tNA\tNA\tNA\tNA"
		myratio = "NA\tNA\tNA\tNA\tNA"
		if myid in dna_abu:
			mydna = dna_abu[myid]
		if myid in rna_abu:
			myrna = rna_abu[myid]
		if myid in ratio_abu:
			myratio = ratio_abu[myid]
		line = re.sub(info[titles["DNA_abundance"]] + "\t" + info[titles["DNA_prevalence"]] + "\t" + info[titles["RNA_abundance"]] + "\t" + info[titles["RNA_prevalence"]] + "\t" + info[titles["RNA-ratio_abundance"]] + "\t" + info[titles["RNA-ratio_prevalence"]], mydna + "\t" + myrna + "\t" + myratio, line)	
		line = re.sub(info[titles["foldChange"]] + "\t" + info[titles["valueChange"]] + "\t" + info[titles["effectSize"]] + "\t" + info[titles["abundance"]] + "\t", "", line)
		open_out.write(line + "\n")
	# foreach file
	open_file.close()
	open_out.close()

# combine_info 


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start combine_abundance_annotation.py -i " + values.i + " ####\n")
	
	### collect abundance info ###
	sys.stderr.write("Get abundance info ......starting\n")
	dna_abu = extract_abundance_info (values.a)
	rna_abu = extract_abundance_info (values.b)
	ratio_abu = extract_abundance_info (values.c)
	combine_info (dna_abu, rna_abu, ratio_abu, values.i, values.o)
	sys.stderr.write("Get abundance info ......done\n")

	sys.stderr.write("### Finish combine_abundance_annotation.py ####\n\n\n")

# end: main

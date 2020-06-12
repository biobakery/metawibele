#!/usr/bin/env python

"""
MetaWIBELE: gene_abundance module
Get gene abundance for each sample

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
import re
import argparse


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Get gene abundance for each sample
"""

def get_args ():
	parser = argparse.ArgumentParser()
	parser.add_argument('-w', help='the working directory', required=True)
	parser.add_argument('-u', help='the input files', required=True)
	parser.add_argument('-r', help='the reference file', required=True)
	parser.add_argument('-t', help='the number of threads', required=True)
	parser.add_argument('-s', help='the sample name', required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# Mapping with bowtie2
#==============================================================
def bowtie2_mapping (work_dir, ref_seq, input_file, thread, sample):
	# Prepare directory
	print(">>Running botwie2...")
	os.system("mkdir -p " + work_dir)
	os.chdir(work_dir)

	# Index reference sequences
	base_name = os.path.basename(ref_seq)
	dir_name = os.path.dirname(ref_seq)
	extension = re.search("\.([^\.]+)$", base_name)
	extension = extension.group(1)
	base_name = re.sub("." + extension + "$", "", base_name)
	ref_index = os.path.join(dir_name, base_name)
	os.system("ln -s " + ref_index + "*" + "  .")

	# Mapping
	sam = sample + ".sam"
	print("bowtie2 -x " + base_name + " -U " + input_file + " --threads " + thread + " -S " + sam + "\n") 
	os.system("bowtie2 -x " + base_name + " -U " + input_file + " --threads " + thread + " -S " + sam) 
	
	return sam
# function bowtie2_mapping


#==============================================================
# Extract reads count info
#==============================================================
def extract_count_info (work_dir, sam_file, ref_seq, thread):
	os.chdir(work_dir)
	
	# filter out low quality mapping reads
	ann_file = re.sub(".fasta", ".saf.gtf", ref_seq)
	ann_file = re.sub(".fna$", ".saf.gtf", ref_seq)
	bam_file = re.sub(".sam", ".bam", sam_file)
	flt_bam = re.sub(".sam", ".flt.bam", sam_file)
	sort_bam = re.sub(".sam", ".sort.bam", sam_file)
	#print("\n>>Filtering low mappping quality...")
	print("samtools view -b " + sam_file + " > " + bam_file)
	#print("samtools view -bq " + q_cutoff + " " + bam_file + " > " + flt_bam)
	print("samtools sort " + bam_file + " -o " + sort_bam)
	print("samtools index " + sort_bam)
	os.system("samtools view -b " + sam_file + " > " + bam_file)
	#os.system("samtools view -bq " + q_cutoff + " " + bam_file + " > " + flt_bam)
	os.system("samtools sort " + bam_file + " -o  " + sort_bam)
	os.system("samtools index " + sort_bam)
	os.system("rm -f " + sam_file)
	os.system("rm -f " + bam_file)

	# Extract abundance
	bed = re.sub(".bam", ".bed", sort_bam)
	print("\n>>Get abundance...")
	print("featureCounts -F SAF -T " + thread + " -a  " + ann_file +  " -o " +  bed + " " + sort_bam)	
	os.system("featureCounts -F SAF -T " + thread + " -a " + ann_file +  " -o " +  bed + " " + sort_bam)	
	#print("bedtools multicov -bams " + sort_bam + " -bed " + bed_file)	
	#os.system("bedtools multicov -bams " + sort_bam + " -bed " + bed_file)	
# extract_count_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args()

	sys.stderr.write("### Start gene_abundance.py -r " + values.r + " ####\n")
	
	### Mapping ###
	sys.stderr.write("Bowtie2 mapping......starting\n")
	sam_file = bowtie2_mapping (values.w, values.r, values.u, values.t, values.s)
	sys.stderr.write("Bowtie2 mapping......done\n")
	
	### Getting abundance ###
	sys.stderr.write("\nExtract gene abundance......starting\n")
	extract_count_info (values.w, sam_file, values.r, values.t)
	sys.stderr.write("\nExtract gene abundance......done\n")
	
	sys.stderr.write("### Finish gene_abundance.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()
	

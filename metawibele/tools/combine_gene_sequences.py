#!/usr/bin/env python
##########################################################################
# Function: Combine genes sequences from different samples
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 07/07/2019
##########################################################################
import sys
import os
import re
import argparse

from metawibele import utilities

def get_args():	
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', help='input the path of sequences and gff files', required=True)
	parser.add_argument('-e', help='input the extentsion of sequence file', default="ffn")
	parser.add_argument('-o', help='output sequence file', required=True)
	values = parser.parse_args()

	return values

#==============================================================
# collect sequences
#==============================================================
def collect_sequence (gene_path, extension, outfile):	
	sampleid = {}
	filelist = utilities.find_files(gene_path, extension, None)
	open_out = open(outfile, "w")
	outfile1 = re.sub(".fna", "_protein_coding.fna", outfile)
	#open_out1 = open(outfile1, "w")
	for myfile in filelist:
		sample = myfile
		mym = re.search("([^\/]+)$", sample)
		sample = mym.group(1)
		sample = re.sub("." + extension, "", sample)
		mygff = re.sub("." + extension, ".gff", myfile)
		
		# collect protein-coding IDs
		gffs = {}
		open_gff = open(mygff, "r")
		for line in open_gff:
			line = line.strip()
			if not len(line):
				continue
			if re.search("^\#", line):
				continue
			if re.search("^\>", line):
				break
			info = line.split("\t")
			if info[2] == "CDS":	# protein-coding genes
				mym = re.search("^ID=([^\;]+)", info[-1])
				gffs[mym.group(1)] = ""
				# debug
				#print("Protein-coding gene:\t" + mym.group(1))
			else:
				# debug
				print("Not CDS\t" + info[2])
		# foreach line
		open_gff.close()
		
		# output sequences
		open_file = open(myfile, "r")
		flag = 0
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				continue
			if re.search("^\>", line):	# sequence id
				mym = re.search("\>([^\_]+)", line)
				myid = mym.group(1)
				sampleid[sample] = myid
				mym = re.search("\>([\S]+)", line)
				mygene = mym.group(1)
				line = re.sub(myid, sample, line)
				open_out.write(line + "\n")
				if mygene in gffs:
					flag = 1
					#open_out1.write(line + "\n")
				else:
					# debug
					print("Not protein coding sequences!\t" + mygene + "\t" + line)
					flag = 0
				continue
			else:
				open_out.write(line + "\n")
				if flag == 1:
					#open_out1.write(line + "\n")
		# foreach line
		open_file.close()
	# foeach sample
	open_out.close()
	#open_out1.close()
	return sampleid
# function collect_sequence


#==============================================================
# output sequence info
#==============================================================
def output_info (sampleid, outfile):
	open_file = open(outfile, "w")
	open_file.write("sample\tID\n")
	for item in sorted(sampleid.keys()):
		open_file.write(item + "\t" + sampleid[item] + "\n")
	#foreachline
	open_file.close()
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start combine_gene_sequences.py -p " + values.p + " ####\n")
	
	### collect sequence info ###
	sys.stderr.write("Get sequence info ......starting\n")
	collect_sequence (values.p, values.e, values.o)
	sys.stderr.write("Get sequence info ......done\n")
	
	sys.stderr.write("### Finish combine_gene_sequences.py ####\n\n\n")

# end: main

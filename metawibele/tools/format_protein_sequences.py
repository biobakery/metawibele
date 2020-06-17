#!/usr/bin/env python

"""
MetaWIBELE: format_protein_sequences module
Combine protein sequences from different samples

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

from metawibele import utilities


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', help='input the path of protein sequence files', required=True)
	parser.add_argument('-q', help='input the path of protein files including partial info', required=True)
	parser.add_argument('-e', help='specify the extension of protein sequence file', default="faa")
	parser.add_argument('-o', help='output sequence file', required=True)
	parser.add_argument('-m', help='output gene info file', required=True)
	values = parser.parse_args()

	return values


# ==============================================================
# collect sequences
# ==============================================================
def collect_sequence (ann_path, extension, partial_path, outfile):
	filelist = utilities.find_files(ann_path, extension, None)
	open_out = open(outfile, "w")
	outfile2 = re.sub(".faa", ".abnormal_seq.faa", outfile)
	outfile2 = re.sub(".fasta", ".abnormal_seq.fasta", outfile)
	open_out2 = open(outfile2, "w")
	gff = {}
	types = {}
	partial = {}
	for myfile in filelist:
		sample = myfile
		mym = re.search("([^\/]+)$", sample)
		sample = mym.group(1)
		sample = re.sub("." + extension, "", sample)

		# collect gff info that is corresponded the sequences file
		contigs = {}
		mapping = {}
		mygff = re.sub("." + extension, ".gff", myfile)
		if not os.path.isfile(mygff):
			print("Gff file doesn't exist!\t" + mygff)
			continue
		open_gff = open(mygff, "r")
		for line in open_gff:
			line = line.strip()
			if not len(line):
				continue
			if re.search("^#", line):
				if re.search("^##sequence-region", line):
					mym = re.search("##sequence-region\s+([\S]+)\s+([\d]+)\s+([\d]+)", line)
					tmp_contig = mym.group(1)
					contigs[tmp_contig] = str(mym.group(2)) + "\t" + str(mym.group(3))
				if re.search("^# Sequence Data", line):
					mytmp = line.split(";")
					mym = re.sub("\"", "", mytmp[-1])
					mym = re.search("^seqhdr\=([\S]+)[\s]+[\s\S]+len\=([\d]+)", mym)
					tmp_contig = mym.group(1)
					contigs[tmp_contig] = str(1) + "\t" + str(mym.group(2))
				continue
			if re.search("^>", line):
				break
			info = line.split("\t")
			feature = info[2]
			start = info[3]
			stop = info[4]
			strand = info[6]
			desc = info[8]
			myinfo = desc.split(";")
			myid = re.search("ID\=([^\;]+)", desc)
			myid = myid.group(1)
			gene_name = "NA"
			gene_id = "NA"
			gene_num = "NA"
			sample_id = "NA"
			contig_id = info[0]
			contig_len = "NA\tNA"
			if contig_id in contigs:
				contig_len = str(contigs[contig_id])
			for item in myinfo:
				if re.search("locus_tag=", item):
					mym = re.search("locus_tag=([\S]+)", item)
					gene_id = mym.group(1)
				if re.search("Name=", item):
					mym = re.search("Name=([\S]+)", item)
					gene_name = mym.group(1)
			# foreach item
			if not re.search("locus_tag=", desc):
				gene_id = sample + "_" + re.sub("_", "-", myid) 
			if not re.search("Name=", desc):
				gene_name = sample + "_" + re.sub("_", "-", myid) 
			if re.search("\_", gene_id):
				mym = re.search("^([^\_]+)\_([\S]+)", gene_id)
				sample_id = mym.group(1)
				gene_num = mym.group(2)
			if not re.search("locus_tag=", desc):
				gene = gene_id
				sample_id = sample
			else:
				gene = sample + "_" + gene_num
			contig = sample + "_contig_" + contig_id
			if feature == "gene":
				if not sample in gff:
					gff[sample] = {}
				if not gene_id in gff[sample]:
					gff[sample][gene_id] = gene + "\t" + gene_id + "\t" + gene_name + "\t" + start + "\t" + stop + "\t" + strand + "\n" + contig + "\t" + contig_id + "\t" + contig_len + "\n" + sample + "\t" + sample_id
			if feature != "gene" and feature != "mRNA":
				if not sample in types:
					types[sample] = {}
				if not gene_id in types[sample]:
					types[sample][gene_id] = feature
				if feature == "CDS":
					new_id = contig_id + "\t" + start + "\t" + stop + "\t" + strand
					mapping[new_id] = gene + "\t" + gene_id
					if not re.search("locus_tag=", desc):
						if not sample in gff:
							gff[sample] = {}
						if not gene_id in gff[sample]:
							gff[sample][gene_id] = gene + "\t" + gene_id + "\t" + gene_name + "\t" + start + "\t" + stop + "\t" + strand + "\n" + contig + "\t" + contig_id + "\t" + contig_len + "\n" + sample + "\t" + sample_id
				# foreach line
		open_gff.close()

		# collect sequences from prodigal results including pratial info
		myfile1 = re.sub(ann_path, partial_path, myfile)
		open_file = open(myfile1, "r")
		AA_seq = {}
		myname = ""
		flag = 0
		hit_num = 0
		total_num = 0
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				continue
			if re.search("^>", line):  # sequence id
				total_num = total_num + 1
				line = re.sub("^>", "", line)
				info = line.split(" # ")
				if len(info) < 4:
					# debug
					print("No info!\t" + myfile1 + "\t" + line)
					continue
				myref = re.sub("_[\d]+$", "", info[0])
				mystart = info[1]
				mystop = info[2]
				mystrand = "+"
				if info[3] == "-1":
					mystrand = "-"
				myid = myref + "\t" + mystart + "\t" + mystop + "\t" + mystrand
				myname = myid
				flag = 0
				if myid in mapping:
					hit_num = hit_num + 1
					# debug
					#print("Mapping:" + myid + "\t" + mapping[myid])
					gene, gene_id = mapping[myid].split("\t")
					myname = ">" + gene
					if not myname in AA_seq:
						AA_seq[myname] = ""
					mym = re.search("partial=([\d]+)", info[-1])
					mypartial = mym.group(1)
					if not sample in partial:
						partial[sample] = {}
					partial[sample][gene_id] = mypartial
					flag = 1
				else:
					# debug
					print("No mapping info!\t" + line)
				continue
			else:
				if flag == 1:
					if myname in AA_seq:
						AA_seq[myname] = AA_seq[myname] + line
		# foreach line
		open_file.close()
	
		if hit_num != total_num:
			open_file = open(myfile, "r")
			myname = ""
			for line in open_file.readlines():
				line = line.strip()
				if not len(line):
					continue
				if re.search("^>", line):  # sequence id	
					mym = re.search("^>([^\_]+)\_([\S]+)", line)
					sample_id = mym.group(1)
					gene_num = mym.group(2)
					gene = sample + "_" + gene_num
					mym = re.search("^>([\S]+)", line)
					gene_id = mym.group(1)
					myname = ">" + gene
					if not myname in AA_seq:
						AA_seq[myname] = ""
					else:
						myname = ""
						continue
					if not sample in partial:
						partial[sample] = {}
					partial[sample][gene_id] = "00"
				else:
					if myname in AA_seq:
						AA_seq[myname] = AA_seq[myname] + line
			# foreach line
			open_file.close()
		

		for myname in sorted(AA_seq.keys()):
			myseq = AA_seq[myname]
			myseq = re.sub("\*$", "", myseq)
			AA_seq[myname] = myseq
			if re.search("\*", myseq): # terminal codon in CDS
				#print("Abnormal CDS\t" + sample + "\t" + myname)
				open_out2.write(myname + "\n" + AA_seq[myname] + "\n")
				continue
			else:
				open_out.write(myname + "\n" + AA_seq[myname] + "\n")
	# foreach sample
	open_out.close()

	return gff, types, partial

# function collect_sequence


# ==============================================================
# output sequence info
# ==============================================================
def output_info (gff, types, partial, outfile):
	open_file = open(outfile, "w")
	#open_file.write("sample\tsample_id\tcontig\tcontig_id\tcontig_start\tcontig_end\tgene\tgene_id\tgene_name\tstart\tend\tstrand\tfeature\tpartial\n")
	open_file.write("GID\tgene\tgene_name\tstart\tend\tstrand\tfeature\tpartial\tCID\tcontig\tcontig_start\tcontig_end\tSID\tsample\n")
	for sample in sorted(gff.keys()):
		for item in sorted(gff[sample].keys()):
			mytype = "NA"
			if sample in types:
				if item in types[sample]:
					mytype = types[sample][item]
			mypartial = "NA"
			if sample in partial:
				if item in partial[sample]:
					mypartial = partial[sample][item]
			if mytype == "NA":
				print("Unknown feature\t" + sample + "\t" + item)
			if mypartial == "NA":
				print("Unknown partial\t" + sample + "\t" + item)
			tmp = gff[sample][item].split("\n")
			open_file.write(tmp[0] + "\t" + mytype + "\t" + mypartial + "\t" + tmp[1] + "\t" + tmp[2] + "\n")
	# foreach line
	open_file.close()


# output_info


# ==============================================================
###########  Main processing ############
# ==============================================================
def main():	
	### get arguments ###

	values = get_args()

	sys.stderr.write("### Start format_protein_sequences.py -p " + values.p + " ####\n")

	### collect sequence info ###
	sys.stderr.write("Get sequence info ......starting\n")
	gff, types, partial = collect_sequence (values.p, values.e, values.q, values.o)
	sys.stderr.write("Get sequence info ......done\n")

	### Output sample info
	sys.stderr.write("\nOutput sequence info ......starting\n")
	output_info (gff, types, partial, values.m)
	sys.stderr.write("Output sequence info ......done\n")

	sys.stderr.write("### Finish format_protein_sequences.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

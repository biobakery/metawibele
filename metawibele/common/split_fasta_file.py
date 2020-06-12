#!/usr/bin/env python


"""
MetaWIBELE: split_family_abundance module
Split the abundance of families

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


#==============================================================
# split files
#==============================================================
def split_fasta_file (seq, split_num, prefix, output, list_file, mylist):
	seqs = {}
	open_file = open(seq, "r")
	myid = ""
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([\S]+)", line)
			myid = ">" + mym.group(1)
			seqs[myid] = ""
		else:
			seqs[myid] = seqs[myid] + line
	# foreach line
	open_file.close()
	
	total_num = len(seqs.keys())
	chunck = int(total_num / int(split_num))
	mynum = 0
	filenum = 0
	start = 0
	out_list = []
	out_list_file = []
	for myid in seqs.keys():
		if mynum > chunck:	# close a split file
			open_out.close()
			mynum = 0 
		mynum = mynum + 1 
		if mynum == 1: # open a new split file
			filenum = filenum + 1 
			myfile = prefix + ".split" + str(filenum) + ".fasta"
			mydir = output + "/" + "split" + str(filenum)
			os.system("mkdir " + mydir)
			myfile = mydir + "/" + myfile
			open_out = open(myfile, "w")
			out_list.append("split" + str(filenum))
			out_list_file.append(myfile)
			open_out.write(myid + "\n" + seqs[myid] + "\n")
		else:
			open_out.write(myid + "\n" + seqs[myid] + "\n")
	# foreach sequence
	open_out.close()

	# ouput file
	open_list = open(mylist, "w")
	for item in out_list:
		open_list.write(item + "\n")
	open_list.close()

	open_list_file = open(list_file, "w")
	for item in out_list_file:
		open_list_file.write(item + "\n")
	open_list_file.close()

# split_fasta_file


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-i', help='input fasta file', required=True)
	parser.add_argument('-n', help='the number of splited files', required=True)
	parser.add_argument('-p', help='prefix of splited file name', required=True)
	parser.add_argument('-w', help='working directory of output files', required=True)
	parser.add_argument('-o', help='output file of splited file names', required=True)
	parser.add_argument('-l', help='output file of splited flags', required=True)
	values=parser.parse_args()

	sys.stderr.write("### Start split_fasta_file.py -i " + values.i + " ####\n")

	split_fasta_file (values.i, values.n, values.p, values.w, values.o, values.l)

	sys.stderr.write("### Finish split_fasta_file.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

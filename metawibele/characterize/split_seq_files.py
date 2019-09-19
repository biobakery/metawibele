#!/usr/bin/env python

##########################################################################
# Function: split large sequence files
# Author: Yancong Zhang (zhangyc201211@gmail.com)
##########################################################################
import sys
import os
import re
import argparse

#==============================================================
# split files
#==============================================================
def split_fasta_file (seq, number, prefix, output, list_file, mylist):
	open_list_file = open(list_file, "w")
	open_list = open(mylist, "w")
	mynum = 0
	filenum = 0
	start = 0
	open_file = open(seq, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue                      
		if re.search("^>", line):
			mym = re.search(">([\S]+)", line)
			line = ">" + mym.group(1)
			if mynum > int(number):
				mynum = 0
			mynum = mynum + 1
			if mynum == 1: # a new splitted file
				filenum = filenum + 1
				#mym = re.search("^([^\.]+)", seq)
				#myname = mym.group(1)
				myfile = prefix + ".split" + str(filenum) + ".fasta"
				mydir = output + "/" + "split" + str(filenum)
				os.system("mkdir " + mydir)
				myfile = mydir + "/" + myfile
				if filenum > (int(start) + 1):
					open_out.close()
				open_list_file.write(myfile + "\n")
				open_list.write("split" + str(filenum) + "\n")
				open_out = open(myfile, "w")
		open_out.write(line + "\n")
	# foreach line
	open_file.close()
	open_list.close()
	open_list_file.close()

# split_fasta_file


#==============================================================
###########  Main processing ############
#==============================================================
if __name__ == '__main__':
	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-i', help='total sequence file', required=True)
	parser.add_argument('-n', help='number of sequences in each splitted file', default=1000)
	parser.add_argument('-p', help='the prefix of splited file name', required=True)
	parser.add_argument('-o', help='the output folder', required=True)
	parser.add_argument('-l', help='the list file of output files', required=True)
	parser.add_argument('-s', help='the list file of split names', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start split_seq_files.py -i " + values.i + " ####\n")
	
	### split files ###
	sys.stderr.write("Split files ......starting\n")
	split_fasta_file (values.i, values.n, values.p, values.o, values.l, values.s)
	sys.stderr.write("Split files ......done\n")

	sys.stderr.write("### Finish split_seq_files.py ####\n\n\n")

# end: main

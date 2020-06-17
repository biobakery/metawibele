#!/usr/bin/env python

"""
MetaWIBELE: psortb_protein module
Extract predicted peptides based on PSORTb results

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


try:
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Extract predicted peptides based on PSORTb results
"""

def get_args ():
	parser = argparse.ArgumentParser()
	parser.add_argument('-e', "--extension",
	                    help='file extension',
	                    required=True,
						default="psortb.gram_positive.out.txt")
	parser.add_argument('-p', "--path",
	                    help='input the path of PSORTb file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect psortb info
#==============================================================
def extract_psortb_info (extension, psortb_path):
	filelist = utilities.find_files(psortb_path, extension, None)
	for myfile in filelist:
		# gram+
		#myfile = psortb_path + "/" + samplelist + "/" + samplelist + ".psortb.gram_positive.out.txt"
		if not os.path.isfile(myfile):
			print ("File not exist!\t" + myfile)
		else:
			print ("OK!\t" + myfile)
			open_file = open(myfile, "r")
			myid = ""
			out_p = []
			flag = 0
			for line in open_file:
				line = line.strip()
				if not len(line):
					continue
				if re.search("SeqID:", line):
					mym = re.search("SeqID:\s+([\S]+)", line)
					myid = mym.group(1)
					continue
				# id
				if re.search("Final", line):
					flag = 1
					continue
				if flag == 1:
					if not len(line):
						continue
					line = re.sub("\s+", "\t", line)
					info = line.split("\t")
					mypredict = "NA"
					myscore = 0
					if not re.search("[\S]+", info[0]):
						mypredict = info[1]
						myscore = info[-1]
					else:
						mypredict = info[0]
						myscore = info[-1]
					if mypredict == myscore:
						myscore = 0
					flag = 0
					if re.search("Unknown", mypredict):
						mypredic = "Unknown"
					myscore = re.sub("\s+", "", str(myscore))
					if re.search("[a-zA-Z]+", myscore):
						myscore = 0
					#open_out.write(myid + "\t" + mypredict + "\t" + str(myscore) + "\n")
					out_p.append(myid + "\t" + mypredict + "\t" + str(myscore))
			# foreach line
			myout = re.sub(".txt", ".location.tsv", myfile)
			open_out = open(myout, "w")
			open_out.write("name\ttype\tscore\n")
			for item in out_p:
				open_out.write(item + "\n")
			open_file.close()
		# if file exist

		# gram-
		#myfile = psortb_path + "/" + samplelist + "/" + samplelist + ".psortb.gram_negtive.out.txt"
		myfile1 = re.sub("psortb.gram_positive.out.txt", "psortb.gram_negative.out.txt", myfile)
		if not os.path.isfile(myfile1):
			print ("File not exist!\t" + myfile1)
		else:
			print ("OK!\t" + myfile1)
			open_file = open(myfile1, "r")
			out_n = []
			myid = ""
			flag = 0	
			for line in open_file:
				line = line.strip()
				if not len(line):
					continue
				if re.search("SeqID:", line):
					mym = re.search("SeqID:\s+([\S]+)", line)
					myid = mym.group(1)
					continue
				# id
				if re.search("Final", line):
					flag = 1
					continue
				if flag == 1:
					if not len(line):
						continue
					line = re.sub("\s+", "\t", line)
					info = line.split("\t")
					mypredict = "NA"
					myscore = 0
					if not re.search("[\S]+", info[0]):
						mypredict = info[1]
						myscore = info[-1]
					else:
						mypredict = info[0]
						myscore = info[-1]
					if mypredict == myscore:
						myscore = 0
					flag = 0
					if re.search("Unknown", mypredict):
						mypredic = "Unknown"
					myscore = re.sub("\s+", "", str(myscore))
					if re.search("[a-zA-Z]+", myscore):
						myscore = 0
					#open_out.write(myid + "\t" + mypredict + "\t" + str(myscore) + "\n")
					out_n.append(myid + "\t" + mypredict + "\t" + str(myscore))
			# foreach line
			#myout = re.sub(".txt", ".location.tsv", myfile)
			myout = re.sub("gram_negative.out.txt", "gram_negative.out.location.tsv", myfile1)
			open_out = open(myout, "w")
			open_out.write("name\ttype\tscore\n")
			for item in out_n:
				open_out.write(item + "\n")
			open_out.close()
		# if file exist

		# archaea
		#myfile = psortb_path + "/" + samplelist + "/" + samplelist + ".psortb.archaea.out.txt"
		myfile1 = re.sub("psortb.gram_positive.out.txt", "psortb.archaea.out.txt", myfile)
		if not os.path.isfile(myfile1):
			print ("File not exist!\t" + myfile1)
			continue
		print ("OK!\t" + myfile1)
		open_file = open(myfile1, "r")
		myid = ""
		out_a = []
		flag = 0
		for line in open_file:
			line = line.strip()
			if not len(line):
				continue
			if re.search("SeqID:", line):
				mym = re.search("SeqID:\s+([\S]+)", line)
				myid = mym.group(1)
				continue
			# id
			if re.search("Final", line):
				flag = 1
				continue
			if flag == 1:
				line = re.sub("\s+", "\t", line)
				info = line.split("\t")
				mypredict = "NA"
				myscore = 0
				if not re.search("[\S]+", info[0]):
					mypredict = info[1]
					myscore = info[-1]
				else:
					mypredict = info[0]
					myscore = info[-1]
				if mypredict == myscore:
					myscore = 0
				flag = 0
				if re.search("Unknown", mypredict):
					mypredic = "Unknown"
				myscore = re.sub("\s+", "", str(myscore))
				if re.search("[a-zA-Z]+", myscore):
					myscore = 0
				#open_out.write(myid + "\t" + mypredict + "\t" + str(myscore) + "\n")
				out_a.append(myid + "\t" + mypredict + "\t" + str(myscore))
		# foreach line
		myout = re.sub(".txt", ".location.tsv", myfile1)
		open_out = open(myout, "w")
		open_out.write("name\ttype\tscore\n")
		for item in out_a:
			open_out.write(item + "\n")
		open_file.close()
	# foreach samplelist
# function extract_psortb_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start psortb_protein.py -p " + values.path + " ####\n")
	
	### collect PSORTb info ###
	sys.stderr.write("Get localization info ......starting\n")
	extract_psortb_info (values.extension, values.path)
	sys.stderr.write("Get localization info ......done\n")

	sys.stderr.write("### Finish psortb_protein.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

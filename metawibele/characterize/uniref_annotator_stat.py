#!/usr/bin/env python

"""
MetaWIBELE: uniref_annotator_stat module
Summary the UniRef annotation

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
	from metawibele import config
	from metawibele import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


#==============================================================
# arguments
#==============================================================
description = """
Function: Summary the UniRef annotation
"""

def get_args():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("-s", "--seq",
						help="Query file used by Diamond",
    					required=True)
	parser.add_argument("-d", "--hit",
						help="Results of Diamond",
    					required=True)
	parser.add_argument("-o", "--output",
						help="The prefix for output file",
    					required=True)
	args = parser.parse_args()
	return args
# get_args


#==============================================================
# collect query id info
#==============================================================
def collect_query (queryfile):
	query = {}
	open_file = open(queryfile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search("^>([\S]+)", line)
			myid = mym.group(1)
			query[myid] = ""
		else:
			continue
	# foreach line
	open_file.close()
	return query
# collect_query


#==============================================================
# collect hits info
#==============================================================
def collect_hits (hitfile):
	hits = {}
	open_file = open(hitfile, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		qseqid = info[0]
		sseqid = info[1]
		if re.search("\|[\S]+", sseqid):
			mym = re.search("^([^\|]+)\|", sseqid)
			sseqid = mym.group(1)
		pident = info[2]
		qcov = abs(int(info[5]) - int(info[4])+1)*1.0/int(info[3])
		scov = abs(int(info[8]) - int(info[7])+1)*1.0/int(info[6])
		mcov = min(qcov, scov)
		#mcov = qcov
		if not qseqid in hits:
			hits[qseqid] = []
		hits[qseqid].append(sseqid + "\t" + str(pident) + "\t" + str(qcov) + "\t" + str(mcov) + "\t" + str(info[-1]))
	# foreach line
	open_file.close()
	return hits
# function collect_hits


#==============================================================
# stat info
#==============================================================
def stat_annotation (query, hits, identity, qcov, mcov, outfile): 
	high_id = {}
	high_cov = {}
	stat_id = {}
	stat_cov = {}
	stat = {}
	identity = float(identity) * 100
	#outfile1 = re.sub(".tsv", ".simple.tsv", outfile)
	open_out = open(outfile, "w")
	#open_out1 = open(outfile1, "w")
	#open_out.write(utilities.PROTEIN_ID + "\tsubject\tidentity\tcoverage\tevalue\ttype\n")
	open_out.write(utilities.PROTEIN_ID + "\tsubject\tidentity\tquery_coverage\tmutual_coverage\tevalue\tquery_type\tmutual_type\n")
	for myid in sorted(query.keys()):
		if not myid in hits:	# no hit
			# debug
			print("No UniRef90 mapping information\t" + myid)
			mystr = myid + "\tNA\tNA\tNA\tNA\tno_hit\tno_hit"
			mystr1 = myid + "\tNA\tNA\tNA\tNA\tNA\tno_hit\tno_hit"
			if not "no_hit" in stat:
				stat["no_hit"] = 1
			else:
				stat["no_hit"] = stat["no_hit"] + 1
			#open_out.write(mystr + "\n")
			open_out.write(mystr1 + "\n")
			continue
		myhit = hits[myid][0]	# default to select the first hit
		myhit_flag = 0
		for item in hits[myid]:
			info = item.split("\t")
			myiden = info[1]
			myqcov = info[2]
			mymcov = info[3]
			flag_id = 0
			flag_cov = 0
			if float(myiden) >= float(identity):
				flag_id = 1
			if float(mymcov) >= float(mcov):
				flag_cov = 1
			if flag_id == 1 and flag_cov == 1:
				myhit = item
				myhit_flag = 1
				break
		# check each hit
		if myhit_flag == 0:
			for item in hits[myid]:
				info = item.split("\t")
				myiden = info[1]
				myqcov = info[2]
				mymcov = info[3]
				flag_id = 0
				flag_cov = 0
				if float(myiden) >= float(identity):
					flag_id = 1
				if float(myqcov) >= float(qcov):
					flag_cov = 1
				if flag_id == 1 and flag_cov == 1:
					myhit = item
					break
			# foreach hit
		# if not found mutual_coverage hit
		info = myhit.split("\t")
		myiden = info[1]
		myqcov = info[2]
		mymcov = info[3]
		mystr = myid + "\t" + info[0] + "\t" + info[1] + "\t" + myqcov + "\t" + info[4]
		mystr1 = myid + "\t" + info[0] + "\t" + info[1] + "\t" + myqcov + "\t" + mymcov + "\t" + info[4]
		flag = 0
		if not myiden in stat_id:
			stat_id[myiden] = {}
		stat_id[myiden][myid] = ""
		if not myqcov in stat_cov:
			stat_cov[myqcov] = {}
		stat_cov[myqcov][myid] = ""
		if float(myiden) >= float(identity):
			flag = flag + 1
			if not myqcov in high_id:
				high_id[myqcov] = {}
			high_id[myqcov][myid] = ""
		if float(myqcov) >= float(qcov):
			flag = flag + 1
			if not myiden in high_cov:
				high_cov[myiden] = {}
			high_cov[myiden][myid] = ""
		if flag == 2:	
			mytype = "high_confidence"
			if not "high_confidence" in stat:
				stat["high_confidence"] = 1
			else:
				stat["high_confidence"] = stat["high_confidence"] + 1
		else:
			mytype = "low_confidence"
			if not "low_confidence" in stat:
				stat["low_confidence"] = 1
			else:
				stat["low_confidence"] = stat["low_confidence"] + 1
		mutual_type = "low_confidence"
		if float(myiden) >= float(identity) and float(mymcov) >= float(mcov):
			mutual_type = "high_confidence"	
		#open_out.write(mystr + "\t" + mytype + "\n")	
		open_out.write(mystr1 + "\t" + mytype + "\t" + mutual_type + "\n")	
	# foreach line
	open_out.close()
	#open_out1.close()

	## output stat
	stat_high_id = re.sub(".tsv", ".high_identity.tsv", outfile)
	stat_high_cov = re.sub(".tsv", ".high_coverage.tsv", outfile)
	out_id = re.sub(".tsv", ".identity.tsv", outfile)
	out_cov = re.sub(".tsv", ".coverage.tsv", outfile)
	stat_all = re.sub(".tsv", ".hits.tsv", outfile)
	
	open_out1 = open(stat_high_id, "w")
	open_out1.write("coverage\tnumber\n")
	for mycov in sorted(high_id.keys(), key=float):
		mynum = len(high_id[mycov].keys())
		open_out1.write(str(mycov) + "\t" + str(mynum) + "\n")
	# foreach coverage
	open_out1.close()
	open_out2 = open(stat_high_cov, "w")
	open_out2.write("identity\tnumber\n")
	for myiden in sorted(high_cov.keys(), key=float):
		mynum = len(high_cov[myiden].keys())
		open_out2.write(str(myiden) + "\t" + str(mynum) + "\n")
	# foreach identity
	open_out2.close()
	
	open_out3 = open(stat_all, "w")
	open_out3.write("type\tnumber\n")
	for mytype in sorted(stat.keys()):
		open_out3.write(str(mytype) + "\t" + str(stat[mytype]) + "\n")
	# foreach type
	open_out3.close()
	
	open_out4 = open(out_id, "w")
	open_out4.write("identity\tnumber\n")
	for myid in sorted(stat_id.keys(), key=float):
		mynum = len(stat_id[myid].keys())
		open_out4.write(str(myid) + "\t" + str(mynum) + "\n")
	# foreach coverage
	open_out4.close()
	open_out5 = open(out_cov, "w")
	open_out5.write("coverage\tnumber\n")
	for mycov in sorted(stat_cov.keys(), key=float):
		mynum = len(stat_cov[mycov].keys())
		open_out5.write(str(mycov) + "\t" + str(mynum) + "\n")
	# foreach identity
	open_out5.close()
# stat_annotation 


#==============================================================
###########  Main processing ############
#==============================================================
def main():

	### get arguments ###
	values = get_args()


	sys.stderr.write("### Start uniref_annotator_stat.py ####\n")
	
	
	### collect info ###
	sys.stderr.write("Get info ......starting\n")
	query = collect_query (values.seq)
	hits = collect_hits (values.hit)
	sys.stderr.write("Get info ......done\n")
	
	### Stat info ###
	sys.stderr.write("Statistic annotation info ......starting\n")
	stat_annotation(query, hits, config.diamond_identity, config.diamond_query_coverage, config.diamond_mutual_coverage, values.output)
	sys.stderr.write("Statistic annotation info ......done\n")


	sys.stderr.write("### Finish uniref_annotator_stat.py ####\n\n")

# end: main

if __name__ == '__main__':
	main()

#!/usr/bin/env python

"""
MetaWIBELE: extract_cluster_CD-hit module
Extract the cluster information from CD-hit results

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
# collect clustering info
#==============================================================
def collect_clustering_info (clust_file):	# combined_genes.sorted.centroid.fna.clstr
	cluster = {}
	length = {}
	if not os.path.isfile(clust_file):
		print ("File not exist!\t" + clust_file)
	else:
		open_file = open(clust_file, "r")
		myclust = ""
		tmp = []
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				continue
			if re.search("^>", line): # new cluster
				if len(tmp) != 0:
					for item in tmp:
						cluster[myclust][item] = ""
				tmp = []
				continue
			mym = re.search('>([^\.]+)\.', line)
			myid = mym.group(1)
			if re.search("\*$", line): # representative
				mym = re.search("\t([\d]+)", line)
				mylen = mym.group(1)
				if not myid in cluster:
					cluster[myid] = {}
				myclust = myid
				length[myid] = mylen
			else: # member
				#if not myclust in cluster:
				#	cluster[myclust] = {}
				#cluster[myclust][myid] = ""
				tmp.append(myid)
		# foreach line
		if len(tmp) != 0:
			for item in tmp:
				cluster[myclust][item] = ""
		open_file.close()
	# exist cluster file
	return cluster, length
# collect_clustering_info 


#==============================================================
# output info
#==============================================================
def output_info (cluster, length, outfile): 
	open_file = open(outfile, "w")
	order = {}
	for myrep in cluster.keys():
		mysize = len(cluster[myrep].keys())
		mysize = mysize + 1
		if not mysize in order:
			order[mysize] = {}
		order[mysize][myrep] = ""
	# foreach rep

	outfile2 = outfile + ".size.tsv"
	open_file2 = open(outfile2, "w")
	open_file2.write("method\tcluster\tsize\n")
	clust_id = 0
	for mysize in sorted(order.keys(), key=int, reverse=True):
		for myrep in order[mysize].keys():
			clust_id = clust_id + 1
			open_file.write(">" + myrep + ";" + "Cluster_" + str(clust_id) + ";length=" + str(length[myrep]) + ";size=" + str(mysize) + ";cluster=" + str(clust_id) + "\n")
			open_file.write(myrep + "\n")
			open_file2.write("CD-hit\t" + "Cluster_" + str(clust_id) + "\t" + str(mysize) + "\n")
			for member in cluster[myrep].keys():
				open_file.write(member + "\n")
		# forech rep
	# foreach size
	open_file.close()
	open_file2.close()
# output_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-c', help='the clust file of CD-hit', required=True)
	parser.add_argument('-o', help='output cluster file', required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start extract_cluster_CD_hit.py -c " + values.c + " ####\n")
	
	### collect clustering info ###
	sys.stderr.write("Get cluster info ......starting\n")
	cluster, length = collect_clustering_info (values.c)
	sys.stderr.write("Get cluster info ......done\n")
	
	### Output summary info
	sys.stderr.write("\nOutput clustering summary info ......starting\n")
	output_info (cluster, length, values.o)
	sys.stderr.write("Output clustering summary info ......done\n")

	sys.stderr.write("### Finish extract_cluster_CD_hit.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

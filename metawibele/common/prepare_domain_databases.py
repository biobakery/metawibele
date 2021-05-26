#!/usr/bin/env python

"""
MetaWIBELE: prepare_domain_databases module
Extract the Pfam name for each Pfam item 

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

try:
	from metawibele.common import utils
	from metawibele import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
            " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Extract the Pfam name for each Pfam item 
"""

def get_args ():
	parser=argparse.ArgumentParser()
	parser.add_argument('-t', "--type", 
						help='specify the Pfam version for downloading, e.g. Pfam33.1', 
						default="Pfam32.0")
	parser.add_argument('-o', "--output",
						help='the path of output folder', 
						required=True)
	values=parser.parse_args()
	return values
# get_args


#==============================================================
# Download the Pfam databases
#==============================================================
def download_dat (output_path, pfam_version):
	dat_file = os.path.join(output_path, "Pfam-A.hmm.dat.gz")
	pfam2go = os.path.join(output_path, "gene_ontology.txt.gz")
	domain = os.path.join(output_path, "INTERACTION.txt.gz")
	pdb_chain_pfam = os.path.join(output_path, "pdb_chain_pfam.tsv.gz")
	pdb_chain_taxonomy = os.path.join(output_path, "pdb_chain_taxonomy.tsv.gz")

	# download pfam
	os.chdir(output_path)
	mypath = os.path.join("ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/", pfam_version)
	download_cmd = "connect " + mypath
	download_cmd = download_cmd + "\n" + "get Pfam-A.hmm.dat.gz"
	download_cmd = download_cmd + "\n " + "get Pfam-A.hmm.gz"
	download_cmd = download_cmd + "\n" + "get Pfam.version.gz"
	download_cmd = download_cmd + "\n" + "connect " + os.path.join(mypath, "database_files")
	download_cmd = download_cmd + "\n" + "get gene_ontology.txt.gz"
	mydownload = os.path.join(output_path, "download_pfam.lftp")
	open_out = open(mydownload, "w")
	open_out.write(download_cmd + "\n")
	open_out.close()
	os.system("lftp -f " + mydownload)

	# download DOMINE
	mypath = "https://manticore.niehs.nih.gov/Domine-2.0/domine-tables-2.0.zip"
	download_cmd = "wget " + mypath
	download_cmd = download_cmd + "; " + "unzip domine-tables-2.0.zip"
	download_cmd = download_cmd + "; " + "gzip INTERACTION.txt"
	os.system(download_cmd)

	# download SIFT data
	mypath = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/"
	download_cmd = "connect " + mypath 
	download_cmd = download_cmd + "\n" + "get pdb_chain_pfam.tsv.gz"
	download_cmd = download_cmd + "\n" + "get pdb_chain_taxonomy.tsv.gz"
	mydownload = os.path.join(output_path, "download_sift.lftp")
	open_out = open(mydownload, "w")
	open_out.write(download_cmd + "\n")
	open_out.close()
	os.system("lftp -f " + mydownload)
	
	return dat_file


#==============================================================
# Extract annotation info
#==============================================================
def extract_annotation_info (datfile, output_path):	# Pfam-A.hmm.dat
	anns = {}
	myid = ""
	myann = ""
	if not os.path.isfile(datfile):
		sys.exit("Error: pfam file doesn't exit! " + datfile)
	outfile = os.path.join(output_path, "pfam_descriptions.txt")
	outfile1 = re.sub(".txt", ".txt.gz", outfile)
	for line in utils.gzip_bzip2_biom_open_readlines (datfile):
		line = line.strip()
		if not len(line):
			continue
		if re.search("\#=GF\s+AC", line):
			mym = re.search("\#=GF\s+AC\s+([^\.]+)", line)
			myid = mym.group(1)
			continue
		if re.search("\#=GF\s+DE", line):
			mym = re.search("\#=GF\s+DE\s+([\S\s]+)", line)
			myann = mym.group(1)
			if not myid in anns:
				anns[myid] = myann
			else:
				anns[myid] = anns[myid] + ";" + myann
			continue
	# foreach line

	open_file = open(outfile, "w")
	open_file.write("Pfam\tdescription\n")
	for mypfam in sorted(anns.keys()):
		open_file.write(mypfam + "\t" + anns[mypfam] + "\n")
	# foreach Pfam
	open_file.close()
	os.system("gzip " + outfile)

# extract_annotation_info 


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	### get arguments ###
	values = get_args ()

	config.logger.info ("### Start prepare_domain_databases step ####")
	
	output = os.path.abspath(values.output)
	if not os.path.isdir(output):
		os.system("mkdir -p " + output)
	
	### Download ###
	config.logger.info ("Download domain info......starting")
	datfile = download_dat (output, values.type)
	config.logger.info ("Download domain......done")
	
	### Extract ###
	config.logger.info ("Extract info......starting")
	extract_annotation_info (datfile, output)
	config.logger.info ("Extract info......done")
	
	config.logger.info ("### Finish prepare_domain_databases step ####")

# end: main

if __name__ == '__main__':
	main()

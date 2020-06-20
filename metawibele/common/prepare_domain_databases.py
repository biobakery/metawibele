#!/usr/bin/env python
##########################################################################
# Function: Extract the Pfam name for each Pfam item 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 10/04/2019
##########################################################################
import sys
import os
import re
import argparse


try:
	from metawibele.common import utils
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


	sys.stderr.write("### Start prepare_domain_databases.py -o " + values.output + " ####\n")
	
	### Download ###
	sys.stderr.write("Download domain info......starting\n")
	datfile = download_dat (values.output, values.type)
	sys.stderr.write("Download domain......done\n")
	
	### Extract ###
	sys.stderr.write("\nExtract info......starting\n")
	extract_annotation_info (datfile, values.output)
	sys.stderr.write("\nExtract info......done\n")
	
	sys.stderr.write("### Finish  prepare_domain_databases.py  ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

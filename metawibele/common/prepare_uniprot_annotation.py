#!/usr/bin/env python
##########################################################################
# Function: Extract the annotation from UniProt data file
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 01/20/2019
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


#==============================================================
# Download the latest version of uniprot dat file
#==============================================================
def download_dat (output_path):
	dat_file = os.path.join(output_path, "uniprot.dat.gz")
	if os.path.isfile(dat_file):
		print("Already exist file and skip this step: " + dat_file)
		return dat_file

	os.chdir(output_path)
	download_cmd = "connect ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
	download_cmd = download_cmd + "\n" + "get reldate.txt"
	download_cmd = download_cmd + "\n" + "get uniprot_sprot.dat.gz"
	download_cmd = download_cmd + "\n" + "get uniprot_trembl.dat.gz"
	mydownload = os.path.join(output_path, "download_uniref.lftp")
	open_out = open(mydownload, "w")
	open_out.write(download_cmd + "\n")
	open_out.close()
	os.system("lftp -f " + mydownload)

	# combine Swiss-Prot and TrEMBLE data
	combine_cmd = "less uniprot_sprot.dat.gz > uniprot.dat"
	combine_cmd = combine_cmd + "; less uniprot_trembl.dat.gz >> uniprot.dat"
	combine_cmd = combine_cmd + "; gzip uniprot.dat"
	os.system(combine_cmd)

	return dat_file


#==============================================================
# Extract annotation info
#==============================================================
def extract_annotation_info (datfile, output_path):
	title = "Entry\tEntry name\tGene names\tProtein names\tOrganism\tNCBI_TaxID\tLength\tGene ontology (biological process)\tGene ontology (molecular function)\tGene ontology (cellular component)\tCross-reference (KEGG)\tCross-reference (KO)\tCOG\tInterpro\tTaxonomic lineage (ALL)\tSubcellular location [CC]\tTransmembrane\tSignal peptide\tCross-reference (Pfam)"
	outfile = os.path.join(output_path, "uniprot_annotation.tsv")
	outfile1 = os.path.join(output_path, "uniprot_annotation.tsv.gz")
	if not os.path.isfile(datfile):
		sys.exit("Error: please download the uniprot dat file!")
	if os.path.isfile(outfile) and not os.path.isfile(outfile1):
		os.system("gzip " + outfile)
		print("Already exist file and skip this step: " + outfile1)
		return outfile1
	if os.path.isfile(outfile1):
		print("Already exist file and skip this step: " + outfile1)
		return outfile1
	open_out = open(outfile, "w")
	open_out.write(title + "\n")
	entry = "NA"
	entry_name = "NA"
	prot_name = "NA"
	gene_name = "NA"
	organism = "NA"
	taxa = "NA"
	length = "NA"
	go_bp = "NA"
	go_mf = "NA"
	go_cc = "NA"
	kegg = "NA"
	ko = "NA"
	cog = "NA"
	interpro = "NA"
	tax = "NA"
	sub = "NA"
	trans = "NA"
	signal = "NA"
	pfam = "NA"
	for line in utils.gzip_bzip2_biom_open_readlines (datfile):
		line = line.strip()
		if not len(line):
			continue
		if line == "//": # end of one entry
			# output info
			if entry != "NA":
				tax = tax + ";" + organism
				open_out.write(entry + "\t" + entry_name + "\t" + gene_name + "\t" + prot_name + "\t" + organism + "\t" + taxa + "\t" + length + "\t" + go_bp + "\t" + go_mf + "\t" + go_cc + "\t" + kegg + "\t" + ko + "\t" + cog + "\t" + interpro+ "\t" + tax + "\t" + sub + "\t" + trans + "\t" + signal + "\t" + pfam + "\n")
			entry = "NA"
			entry_name = "NA"
			prot_name = "NA"
			gene_name = "NA"
			organism = "NA"
			taxa = "NA"
			length = "NA"
			go_bp = "NA"
			go_mf = "NA"
			go_cc = "NA"
			kegg = "NA"
			ko = "NA"
			cog = "NA"
			interpro = "NA"
			tax = "NA"
			sub = "NA"
			trans = "NA"
			signal = "NA"
			pfam = "NA"
			continue
		if re.search("^ID\s+", line): # entry name
			mym = re.search("^ID\s+([\S]+)\s+[\s\S]+\s+([\d]+)\s+AA", line)
			entry_name = mym.group(1)
			length = mym.group(2)
			continue
		if re.search("^AC\s+", line): # entry
			mym = re.search("AC\s+([\S]+)", line)
			entry = mym.group(1)
			entry = re.sub(";$", "", entry)
			continue
		if re.search("^GN\s+", line): # entry
			mym = re.search("GN\s+([\S]+[\s\S]+)", line)
			gene_tmp = mym.group(1)
			gene_tmp = re.sub("\s+", "", gene_tmp)
			tmps = gene_tmp.split(";")
			gene_name_tmp = ""
			for tmp in tmps:
				if re.search("\=([^\=]+)", tmp):
					mym = re.search("\=([^\=]+)", tmp)
					gene_name_tmp = gene_name_tmp + mym.group(1) + ";"
			gene_name_tmp = re.sub(";$", "", gene_name_tmp)
			if gene_name == "NA":
				gene_name = gene_name_tmp 
			else:
				gene_name = gene_name + ";" + gene_name_tmp
			continue
		if re.search("^DE\s+[\s\S]+Full\=([\S\s]+)", line): # protein name
			mym = re.search("^DE\s+[\s\S]+Full\=([\S\s]+)", line)
			mydec = mym.group(1)
			mydec = re.sub(";$", "", mydec)
			if prot_name == "NA":
				prot_name = mydec
			else:
				prot_name = prot_name + ";" + mydec
			continue
		if re.search("^OS\s+", line): # organism
			mym = re.search("^OS\s+([\S\s]+)", line)
			mydec = mym.group(1)
			mydec = re.sub("\.$", "", mydec)
			if organism == "NA":
				organism = mydec
			else:
				organism = organism + ";" + mydec
			continue
		if re.search("^OX\s+", line): # NCBI taxonomy ID
			mym = re.search("^OX\s+NCBI_TaxID=([\d]+)", line)
			taxa_info = mym.group(1)
			taxa_info = re.sub(";$", "", taxa_info)
			if taxa == "NA":
				taxa = taxa_info
			else:
				taxa = taxa + ";" + taxa_info
			continue
		if re.search("^OC\s+", line): # taxonomy lineage
			mym = re.search("^OC\s+([\S\s]+)", line)
			mydec = mym.group(1)
			mydec = re.sub("\.$", "", mydec)
			mydec = re.sub(";$", "", mydec)
			if tax == "NA":
				tax = mydec
			else:
				tax = tax + ";" + mydec
			continue
		if re.search("^DR\s+KEGG;", line): # KEGG annotation
			mym = re.search("^DR\s+(KEGG;[\S\s]+)", line)
			mydec = mym.group(1)
			tmp = mydec.split("; ")
			kegg_tmp = tmp[1]
			if kegg == "NA":
				kegg = kegg_tmp
			else:
				kegg = kegg + ";" + kegg_tmp
			continue
		if re.search("^DR\s+KO;", line): # KO annotation
			mym = re.search("^DR\s+(KO;[\S\s]+)", line)
			mydec = mym.group(1)
			tmp = mydec.split("; ")
			ko_tmp = tmp[1]
			if ko == "NA":
				ko = ko_tmp
			else:
				ko = ko + ";" + ko_tmp
			continue
		if re.search("^DR\s+eggNOG;", line): # COG annotation
			mym = re.search("^DR\s+(eggNOG;[\S\s]+)", line)
			mydec = mym.group(1)
			tmp = mydec.split("; ")
			cog_tmp = tmp[1]
			if cog == "NA":
				cog = cog_tmp
			else:
				cog = cog + ";" + cog_tmp
			continue
		if re.search("^DR\s+InterPro;", line): # InterPro annotation
			mym = re.search("^DR\s+(InterPro;[\S\s]+)", line)
			mydec = mym.group(1)
			tmp = mydec.split("; ")
			inter_tmp = tmp[1]
			if interpro == "NA":
				interpro = inter_tmp
			else:
				interpro = interpro + ";" + inter_tmp
			continue
		if re.search("^DR\s+GO;", line): # GO annotation
			mym = re.search("^DR\s+(GO;[\S\s]+)", line)
			mydec = mym.group(1)
			tmp = mydec.split("; ")
			go = tmp[1]
			go_info = tmp[2]
			if re.search("^P:", go_info):
				go_info = re.sub("^P:", "", go_info)
				if go_bp == "NA":
					go_bp = go_info + " [" + go + "]"
				else:
					go_bp = go_bp + ";" + go_info + " [" + go + "]"
			if re.search("^F:", go_info):
				go_info = re.sub("^F:", "", go_info)
				if go_mf == "NA":
					go_mf = go_info + " [" + go + "]"
				else:
					go_mf = go_mf + ";" + go_info + " [" + go + "]"
			if re.search("^C:", go_info):
				go_info = re.sub("^C:", "", go_info)
				if go_cc == "NA":
					go_cc = go_info + " [" + go + "]"
				else:
					go_cc = go_cc + ";" + go_info + " [" + go + "]"
			continue
		if re.search("^DR\s+Pfam;", line): # GO annotation
			mym = re.search("^DR\s+(Pfam;[\S\s]+)", line)
			mydec = mym.group(1)
			tmp = mydec.split("; ")
			pfam_info = tmp[1]
			if pfam == "NA":
				pfam = pfam_info
			else:
				pfam = pfam + ";" + pfam_info
			continue
		if re.search("^CC\s+", line): # CC info
			if re.search("^CC\s+\-\!\-\s+SUBCELLULAR\s+LOCATION", line): # Subcellular
				mym = re.search("^CC\s+\-\!\-\s+SUBCELLULAR\s+LOCATION\s*([\s\S]+)", line)
				myflag = 1
				mydec = mym.group(1)
				if sub == "NA":
					sub = "SUBCELLULAR LOCATION" + mydec
				else:
					sub = sub + ";" + "SUBCELLULAR LOCATION" + mydec
			elif re.search("^CC\s+\-\!\-", line):	# not subcellular
				myflag = 0
			else:
				if myflag == 1:
					mym = re.search("^CC\s+(\S+[\s\S]+)", line)
					mydec = mym.group(1)
					mydec = re.sub("\.$", "", mydec)
					sub = sub + " " + mydec
			continue
		if re.search("^FT\s+", line): # FT info
			if re.search("^FT\s+TRANSMEM\s+", line): # transmembrane
				transflag = 1
				signalflag = 0
				if trans == "NA":
					trans = "TRANSMEM"
				else:
					trans = trans + ";" + "TRANSMEM"
			if re.search("^FT\s+SIGNAL\s+", line):	# signal
				signalflag = 1
				transflag = 0
				if signal == "NA":
					signal = "SIGNAL"
				else:
					signal = signal + ";" + "SIGNAL"
			continue
		# FT info
	# foreach line
	# last entry
	tax = tax + ";" + organism
	open_out.write(entry + "\t" + entry_name + "\t" + prot_name + "\t" + organism + "\t" + taxa + "\t" + length + "\t" + go_bp + "\t" + go_mf + "\t" + go_cc + "\t" + kegg + "\t" + ko + "\t" + interpro + "\t" + cog + "\t" + tax + "\t" + sub + "\t" + trans + "\t" + signal + "\t" + pfam + "\n")
	open_out.close()

	os.system("gzip " + outfile)
# extract_annotation_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():

	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-o', "--output", 
						help='output folder for uniprot annotation info', 
						required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start prepare_uniprot_annotation.py -o " + values.output + " ####\n\n")
	
	### Extraction ###
	sys.stderr.write("Download uniprot info......starting\n")
	datfile = download_dat (values.output)
	sys.stderr.write("Download uniprot info......done\n")
	sys.stderr.write("Extract uniprot info......starting\n")
	extract_annotation_info (datfile, values.output)
	sys.stderr.write("Extract uniprot info......done\n\n")
	
	sys.stderr.write("### Finish  prepare_uniprot_annotation.py  ####\n\n")

# end: main

if __name__ == '__main__':
	main()

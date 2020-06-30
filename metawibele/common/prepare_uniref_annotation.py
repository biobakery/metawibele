#!/usr/bin/env python
##########################################################################
# Function: Extract the annotation info from UniRef fasta file and UniProt annotation file
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 01/20/2019
##########################################################################
import sys
import os
import re
import argparse

try:
	from metawibele import config
	from metawibele.common import utils
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
                 " Please check your install.")

#==============================================================
# Download the latest version of uniref fasta file
#==============================================================
def download_seq (output_path, uniref_type):
	comp_seq = uniref_type + ".fasta.gz"
	uncomp_seq = uniref_type + ".fasta"
	release_note = uniref_type + ".release_note"
	uncomp_file = os.path.join(output_path, os.path.basename(uncomp_seq))
	comp_file = os.path.join(output_path, os.path.basename(comp_seq))
	if os.path.isfile(uncomp_file):
		print("File exist and skipe this step: " + uncomp_file)
		return uncomp_file
	if os.path.isfile(comp_file) and not os.path.isfile(uncomp_file):
		os.sytem("gunzip " + comp_file)
		return uncomp_file

	os.chdir(output_path)
	download_cmd = "connect ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/" 
	download_cmd = download_cmd + "\n" + "get " + comp_seq
	download_cmd = download_cmd + "\n" + "get " + release_note
	mydownload = os.path.join(output_path, "download_uniref.lftp")
	open_out = open(mydownload, "w")
	open_out.write(download_cmd + "\n")
	open_out.close()
	os.system("lftp -f " + mydownload)

    # uncompress fasta file
	combine_cmd = "gunzip " + comp_seq
	os.system(combine_cmd)

	return uncomp_file


#==============================================================
# Extract taxonomy mapping info
#==============================================================
def extract_mapping_info (mapfile):	# uniprot_taxonomy.map.tsv
	maps = {}
	titles = {}
	for line in utils.gzip_bzip2_biom_open_readlines (mapfile):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^Taxon", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		taxa_id = info[titles["Taxon"]]
		taxa_name = info[titles["Scientific_name"]]
		maps[taxa_id] = taxa_name	# TaxID <-> Tax
	# foreach line

	return maps
# extract_mappig_info


#==============================================================
# Extract annotation info
#==============================================================
def extract_annotation_info (output_path, maps):
	anns = {}
	human_pfams = {}
	titles = {}
	datfile = os.path.join(output_path, "uniprot_annotation.tsv.gz")
	if not os.path.isfile(datfile):
		sys.exit("Error: uniprot annotation file (uniprot_annotation.tsv.gz) doesn't exit, please prepare for it")

	items = ["UniProtKB", "Entry_name", "Organism", "Gene_names", "Length", "GO_BP", "GO_MF", "GO_CC", "KEGG", "KEGG-KO", "eggNOG", "Interpro", "Taxonomic_lineage", "Subcellular_location", "Transmembrane", "Signal_peptide", "Pfam"]
	for line in utils.gzip_bzip2_biom_open_readlines (datfile):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^#", line):
			continue
		if re.search("^Entry", line):
			myindex = 0
			while myindex < len(info):
				item = info[myindex]
				if item == "Entry":
					item = "UniProtKB"
				if item == "Entry name":
					item = "Entry_name"
				if item == "Gene names":
					item = "Gene_names"
				if item == "Gene ontology (biological process)":
					item = "GO_BP"
				if item == "Gene ontology (molecular function)":
					item = "GO_MF"
				if item == "Gene ontology (cellular component)":
					item = "GO_CC"
				if re.search("KEGG", item):
					item = "KEGG"
				if re.search("\(KO\)", item):
					item = "KEGG-KO"
				if re.search("Taxonomic", item):
					item = "Taxonomic_lineage"
				if re.search("Subcellular location", item):
					item = "Subcellular_location"
				if re.search("Signal peptide", item):
					item = "Signal_peptide"
				if re.search("Pfam", item):
					item = "Pfam"
				if item == "Protein names":
					item = "Protein_names"
				if item in items:
					titles[myindex] = item
				if "NCBI_TaxID" == item:
					titles[myindex] = item
				if "Protein_names" == item:
					titles[myindex] = item
				myindex = myindex + 1
			# foreach item
			continue
		
		myindex = 0
		mystr = {}
		tax = "NA"
		taxID = "NA" # taxa
		uniprot_id = "NA"
		org = "NA"
		gene = "NA"
		protein = "NA"
		pfam_item = "NA"
		entry_name = info[1]
		while myindex < len(info):
			item = info[myindex]
			if myindex in titles:
				myid = titles[myindex] 
				if item == "":
					item = "NA"
				if myid == "NCBI_TaxID":
					taxID = item
					if taxID in maps:
						tax = maps[taxID]
					myindex = myindex + 1
					continue
				# UniProtKB
				if myid == "UniProtKB" and item != "NA" and item != "":
					uniprot_id = item
				# Gene name
				if myid == "Gene_names" and item != "NA" and item != "":
					gene = item
					gene = re.sub(";$", "", gene)
					gene = re.sub("\{[^\{]+\}", "", gene)
					gene = re.sub("\s+", ";", gene)
				# Protein name
				if myid == "Protein_names" and item != "NA" and item != "":
					protein = item
					protein = re.sub(";$", "", protein)
					protein = re.sub("\{[^\{]+\}", "", protein)
					protein = re.sub("\s+$", "", protein)
				# Organism
				if myid == "Organism" and item != "NA" and item != "":
					org = item
				# Pfam
				if myid == "Pfam" and item != "NA" and item != "":
					pfam_item = item
				mystr[myid] = item
			myindex = myindex + 1
		# foreach item
		mystr_out = ""
		for item in items:
			if item in mystr:
				if mystr_out == "":
					mystr_out = mystr[item]
				else:
					mystr_out = mystr_out + "\t" + mystr[item]
			else:
				if mystr_out == "":
					mystr_out = "NA"
				else:
					mystr_out = mystr_out + "\tNA"
		anns[info[1]] = tax + "\t" + taxID + "\t" + mystr_out
		
		if taxID == "9606":	# human proteins
			if pfam_item != "NA":
				pfam_item = re.sub("\s+", "", pfam_item)
				tmps = pfam_item.split(";")
				for i in tmps: 
					if not i in human_pfams:
						human_pfams[i] = {}
					if uniprot_id != "NA":
						human_pfams[i][uniprot_id] = protein + "\t" + gene
	
	# foreach line
	
	# report human pfams
	outfile = os.path.join(output_path, "uniprot_human_pfam.tsv")
	outfile1 = re.sub(".tsv", ".tsv.gz", outfile)
	open_out = open(outfile, "w")
	open_out.write("Pfam\tOrganism\tNCBI_TaxID\tUniProtKB\tProtein_names\tGene_names\n")
	for myid in sorted(human_pfams.keys()):
		mystr = myid + "\tHomo sapiens (Human)\t9606"
		#mypfam = ";".join(sorted(human_pfams[myid].keys()))
		mystr1 = ""
		mystr2 = ""
		mystr3 = ""
		for i in sorted(human_pfams[myid].keys()):
			mystr1 = mystr1 + i + ";"
			x, y = human_pfams[myid][i].split("\t")
			if x != "NA" and x != "":
				mystr2 = mystr2 + x + ";"
			if y != "NA" and y != "":
				mystr3 = mystr3 + y + ";"
		mystr1 = re.sub(";$", "", mystr1)
		mystr2 = re.sub(";$", "", mystr2)
		mystr3 = re.sub(";$", "", mystr3)
		if mystr1 == "":
			mystr1 = "NA"
		if mystr2 == "":
			mystr2 = "NA"
		if mystr3 == "":
			mystr3 = "NA"
		open_out.write(mystr + "\t" + mystr1 + "\t" + mystr2 + "\t" + mystr3 + "\n")
	open_out.close()
	os.system("gzip " + outfile)

	header = "Rep_Tax\tRep_TaxID\t" + "\t".join(items)
	
	return anns, header
# extract_annotation_info


#==============================================================
# Extract taxonomy info
#==============================================================
def extract_info (ann_info, ann_title, seqfile, uniref_type, output_path):
	title = "ID\tTax\tTaxID\tProtein_names\t" + ann_title 
	outfile = os.path.join(output_path, uniref_type + "_annotation.tsv")
	outfile1 = re.sub(".tsv", ".tsv.gz", outfile)
	if not os.path.isfile(seqfile):
		sys.exit("Error: " + seqfile + " doesn't exit!")
	if os.path.isfile(outfile1):
		print("Already exist file and skip this step: " + outfile1)
		return(outfile1)
	if os.path.isfile(outfile) and not os.path.isfile(outfile1):
		os.system("gzip " + outfile)
		return(outfile1)

	open_file = open(seqfile, "r")
	open_out = open(outfile, "w")
	open_out.write(title + "\n")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			if re.search(">([\S]+)([\S\s]+)n=([\d]+)\s*Tax=([\S\s]+)TaxID=([\d]+)\s*RepID=([\S]+)", line):
				mym = re.search(">([\S]+)([\S\s]+)n=([\d]+)\s*Tax=([\S\s]+)TaxID=([\d]+)\s*RepID=([\S]+)", line)
				myid = mym.group(1)
				mydec = mym.group(2)
				mydec = mydec.strip()
				mynum = mym.group(3)
				mytax = mym.group(4)
				mytax = mytax.strip()
				mytaxID = mym.group(5)
				myrep = mym.group(6)
				myann = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" 
				if myrep in ann_info:
					# debug
					#print("Hit\t" + myrep)
					myann = ann_info[myrep]
				#ann[myid] = myid + "\t" + mytax + "\t" + mytaxID + "\t" + mydec + "\t" + myann
				tmp = myann.split("\t")
				if re.search("_UPI", myid) and tmp[0] == "NA":
					tmp[0] = mytax
					tmp[1] = mytaxID
					myann = "\t".join(tmp)
				mystr = myid + "\t" + mytax + "\t" + mytaxID + "\t" + mydec + "\t" + myann
				open_out.write(mystr + "\n")
			else:
				if re.search(">([\S]+)([\S\s]+)n=([\d]+)\s*Tax=unknown[\s\S]*RepID=([\S]+)", line):
					mym = re.search(">([\S]+)([\S\s]+)n=([\d]+)\s*Tax=(unknown)[\s\S]*RepID=([\S]+)", line)
					myid = mym.group(1)
					mydec = mym.group(2)
					mydec = mydec.strip()
					mynum = mym.group(3)
					mytax = mym.group(4)
					mytax = mytax.strip()
					mytaxID = "NA"
					myrep = mym.group(5)
					myann = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" 
					if myrep in ann_info:
						myann = ann_info[myrep]
					tmp = myann.split("\t")
					if re.search("_UPI", myid) and tmp[0] == "NA":
						tmp[0] = mytax
						tmp[1] = mytaxID
						myann = "\t".join(tmp)
					mystr = myid + "\t" + mytax + "\t" + mytaxID + "\t" + mydec + "\t" + myann
					open_out.write(mystr + "\n")
				else:
					print("Exception\t" + line)
		else:
			continue
	# foreach line
	open_file.close()

	os.system("gzip " + outfile)

	return outfile1
# function extract_info


def report_each_map (uniref_ann, map_type, uniref_type, output_path):
	outfile = os.path.join(output_path, "map_" + map_type + "_" + uniref_type + ".txt")
	outfile1 = re.sub(".txt", ".txt.gz", outfile)
	outs = {}
	titles = {}
	titles1 = {}
	ids = {}
	for line in utils.gzip_bzip2_biom_open_readlines (uniref_ann):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^ID", line) or re.search("^UniRef\t", line) or re.search("^UniRefID\t", line):
			for item in info:
				titles[item] = info.index(item)		
				titles1[info.index(item)] = item
			continue
		info = line.split("\t")
		if "UniRef" in titles:
			myuniref = info[titles["UniRef"]]
		elif "ID" in titles:
			myuniref = info[titles["ID"]]
		myindex = 0 
		while myindex < len(info):
			mykey = titles1[myindex]
			mykey = re.sub("\(", "_", mykey)
			mykey = re.sub("\)", "", mykey)
			myvalue = info[myindex]
			if myvalue == "NA" or myvalue == "na" or myvalue == "NaN":
				myindex = myindex + 1
				continue
			if mykey == map_type:
				if not myvalue in outs:
					outs[myvalue] = {}
				outs[myvalue][myuniref] = ""
			myindex = myindex + 1
    # foreach uniref

	open_out = open(outfile, "w")
	for myid in sorted(outs.keys()):
		myinfo = "\t".join(sorted(outs[myid].keys()))
		open_out.write(myid + "\t" + myinfo + "\n")
	open_out.close()
	
	os.system("gzip " + outfile)


def extract_uniref_maps (uniref_ann, uniref_type, output_path):
	types = ["Protein_names", "Gene_names", "UniProtKB", "Tax", "TaxID", "Rep_Tax", "Rep_TaxID", "GO", "GO_BP", "GO_CC", "GO_MF", "KEGG-KO", "eggNOG", "Pfam"] 
	for mytype in types:
		report_each_map (uniref_ann, mytype, uniref_type, output_path)



#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	parser=argparse.ArgumentParser()
	parser.add_argument('-t', "--type",
						help='specify the type of uniref', 
						choices=['uniref90', 'uniref50'], 
						default='uniref90')
	parser.add_argument('-o', "--output", 
						help='output folder', 
						required=True)
	values=parser.parse_args()


	sys.stderr.write("### Start prepare_uniref_annotation.py -o " + values.output + " ####\n")
	
	### Extraction ###
	sys.stderr.write("Download uniref sequences......starting\n")
	seqfile = download_seq (values.output, values.type)
	sys.stderr.write("Download uniref sequences......done\n")
	maps = extract_mapping_info (config.taxonomy_database)
	sys.stderr.write("Collect uniprot annotations......starting\n")
	ann_info, ann_title = extract_annotation_info (values.output, maps)
	sys.stderr.write("Collect uniprot annotations......done\n")
	
	### Output ###
	sys.stderr.write("\nReport uniref combined info......starting\n")
	uniref_ann = extract_info (ann_info, ann_title, seqfile, values.type, values.output)
	anns_info = {}
	sys.stderr.write("\nReport uniref combined info......done\n")
	sys.stderr.write("\nReport uniref mapping info......starting\n")
	extract_uniref_maps (uniref_ann, values.type, values.output)
	sys.stderr.write("\nReport uniref mapping info......done\n")
	
	sys.stderr.write("### Finish prepare_uniref_annotation.py  ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

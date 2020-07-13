#!/usr/bin/env python

"""
MeteWIBELE: finalize characterization module
1) combine all functional and taxonomic annotations
2) format outputs

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

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Finalize the annotation results 
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--annotation",
	                    help='input pre-summarized annotation file', required=True)
	parser.add_argument('-t', "--taxonomy",
	                    help='input taxonomic annotation file', required=True)
	parser.add_argument('-u', "--uniref",
	                    help='input uniref90 taxonomic annotation file', required=True)
	parser.add_argument('-l', "--list",
	                    help='input list file of annotation files', required=True)
	parser.add_argument('-s', "--source",
	                    help='specify source of annotations',
	                    choices=["protein", "protein_family"],
	                    required=True,
	                    default="protein_family")
	parser.add_argument('-o', "--output",
	                    help='output annotation file',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_cluster_info (clust_file, source):	# discovery_cohort.peptides.clust
	cluster = {}
	open_file = open(clust_file, "r")
	myclust_id = "NA"
	myclust = "NA"
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			mym = re.search("length=([\d]+)", line)
			mylen = mym.group(1)
			mym = re.search("size=([\d]+)", line)
			mysize = mym.group(1)
			mym = re.search("cluster=([\d]+)", line)
			myclust_id = "Cluster_" + mym.group(1)
			if source == "protein_family":
				cluster[myclust_id] = myclust_id + "\t" + "protein_family" + "\t" + "repID=" + myclust + "\n" + "rep_length=" + mylen + "\n" + "cluster_size=" + mysize
			continue
		mym = re.search("([\S]+)", line)
		myid = mym.group(1)
		if source == "protein_family":
			continue
		cluster[myid] = myclust_id + "\t" + "protein_family" + "\t" + "repID=" + myclust + "\n" + "rep_length=" + mylen + "\n" + "cluster_size=" + mysize
	# foreach line
	open_file.close()
	return cluster
# function collect_cluster_info


#==============================================================
# collect taxonomic annotation for protein families
#==============================================================
def collect_taxonomy_annotation (taxa_ann):
	taxonomy = {}
	mapping = {}
	titles = {}
	open_file = open(taxa_ann, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	id_flag = info[0]
	for item in info:
		titles[item] = info.index(item)
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		map_type = info[titles["map_type"]]
		query_type = info[titles["query_type"]]
		mutual_type = info[titles["mutual_type"]]
		identity = info[titles["identity"]]
		query_cov = info[titles["query_coverage"]]
		mutual_cov = info[titles["mutual_coverage"]]
		mydesc = info[titles["detail"]]
		tax = info[titles["Tax"]]
		taxID = info[titles["TaxID"]]
		reptax = info[titles["Rep_Tax"]]
		reptaxID = info[titles["Rep_TaxID"]]
		myname = info[titles["Tax"]]
		msp_name = "NA"
		msp_taxa_name = "NA"
		msp_taxa_id = "NA"
		if myname == "NA" or myname == "root" or myname == "Bacteria" or myname == "unclassified sequences" or myname == "cellular organisms":
			myname = "Unclassified"
		if "MSP_Tax" in titles:
			tax = info[titles["MSP_Tax"]]
		if "MSP_TaxID" in titles:
			taxID = info[titles["MSP_TaxID"]]
		if "MSP_Rep_Tax" in titles:
			reptax = info[titles["MSP_Rep_Tax"]]
		if "MSP_Rep_TaxID" in titles:
			reptaxID = info[titles["MSP_Rep_TaxID"]]
		if "msp_name" in titles:
			msp_name = info[titles["msp_name"]]
		if "msp_taxa_name" in titles:
			msp_taxa_name = info[titles["msp_taxa_name"]]
		if "msp_taxa_id" in titles:
			msp_taxa_id = info[titles["msp_taxa_id"]]
		#org = info[titles["organism"]]
		unirefID = info[titles["unirefID"]]
		uniprot = info[titles["UniProtKB"]]
		taxa_id = info[titles["taxa_id"]]
		taxa_name = info[titles["taxa_name"]]
		taxa_rank = info[titles["taxa_rank"]]
		taxa_lineage = info[titles["taxa_lineage"]]
		if taxa_rank == "Unclassified":
			myname = "Unclassified"
		if map_type == "UniRef90_uncharacterized" or map_type == "UniRef90_characterized":
			map_type = "strong_homology"
		if map_type == "UniRef90_weak_homology":
			map_type = "weak_homology"
		if map_type == "UniRef90_worse_homology":
			map_type = "worse_homology"
		if map_type == "UniRef90_none_homology":
			map_type = "none_homology"
			taxonomy[myid] = "NA\tNA\tNA" 
			mapping[myid] = "NA\t" + map_type + "\tNA"
		else:
			mapping[myid] = unirefID + "\t" + map_type + "\t" + "UniProtKB=" + uniprot + "\n" + "Protein_names=" + mydesc +  "\n" + "query_cov_type=" + query_type + "\n" + "mutual_cov_type=" + mutual_type + "\n" + "identity=" + identity + "\n" + "query_coverage=" + query_cov + "\n" + "mutual_coverage=" + mutual_cov + "\n" + "taxa_id=" + info[titles["TaxID"]] + "\ntaxa_name=" + myname

			taxonomy[myid] = taxa_name + "\t" + taxa_rank + "\t" + "taxa_id=" + taxa_id + "\n" + "taxa_lineage=" + taxa_lineage + "\n" + "LCA_Tax=" + tax + "\n" + "LCA_TaxID=" + taxID + "\n" + "Rep_Tax=" + reptax + "\n" + "Rep_TaxID=" + reptaxID + "\n" + "msp_name=" + msp_name + "\n" + "msp_taxa_name=" + msp_taxa_name + "\n" + "msp_taxa_id=" + msp_taxa_id
	# foreach line
	open_file.close()
	return taxonomy, mapping
# collect_taxonomy_annotation


#==============================================================
# collect functional annotation for protein families
#==============================================================
def collect_annotation (list_file, source):
	annotation = {}
	open_list = open(list_file, "r")
	for myfile in open_list.readlines():
		myfile = myfile.strip()
		if not len(myfile):
			continue
		if re.search("^#", myfile):
			continue
		if not os.path.isfile(myfile):
			print("File not exist!\t" + myfile)
			continue
		open_file = open(myfile, "r")
		mym = re.search(config.basename + "_([\S]+)_proteinfamilies", os.path.basename(myfile))
		method = mym.group(1)
		titles = {}
		titles_item = {}
		if source == "protein_family":
			myflag = utilities.PROTEIN_FAMILY_ID
		else:
			myflag = utilities.PROTEIN_ID
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			if re.search("^" + myflag, line):
				myindex = 0
				while myindex < len(info):
					item = info[myindex]
					titles[item] = myindex 
					titles_item[myindex] = item
					myindex = myindex + 1
				continue
			myid = info[titles[myflag]]
			mytype = info[titles["type"]]
			detail = info[titles["detail"]]

			# collect attributes
			mystr = "NA"
			myindex = 3
			while myindex < len(info):
				if not myindex in titles_item:
					print("No title item info!\t" + myfile + "\t" + myid + "\t" + str(myindex))
					continue
				myname = titles_item[myindex]
				if mystr == "NA":
					mystr = myname + "=" + info[myindex]
				else:
					mystr = mystr + "\n" + myname + "=" + info[myindex]
				myindex = myindex + 1
			mykey = myid + "\t" + method + "\t" + mytype
			if not mykey in annotation:
				annotation[mykey] = []
			annotation[mykey].append(detail + "\t" + mytype + "\t" + mystr)
		# foreach line
		open_file.close()
	# foreach annotation file
	open_list.close()
	return annotation
# collect_annotation


#==============================================================
# finalize annotation 
#==============================================================
def output_attributes (line, mynum, output_fh):	
	tmp = line.split("\n\n")
	myattr = tmp[0]
	tmp2 = tmp[1].split("\n")
	for item in tmp2:
		if not re.search("([^=]+)=([^=]+)", item):
			# debug
			print("Item with format error!\t" + myattr + "\t" + item)
			continue
		mym = re.search("([^=]+)=([^=]+)", item)
		mykey = mym.group(1)
		myvalue = mym.group(2)
		mynum = mynum + 1
		output_fh.write(str(mynum) + "\t" + myattr + "\t" + mykey + "\t" + myvalue + "\n")
	# foreach attribute
	
	return mynum

def finalize_annotation (source, pep_cluster, annotation, taxonomy, mapping, annfile, outfile):
	titles = {}
	flag_type = {}
	flag_clust = {}
	open_in = open(annfile, "r")
	open_out = open(outfile, "w")
	if source == "protein_family":
		open_out.write(utilities.PROTEIN_FAMILY_ID + "\tannotation\tfeature\tcategory\tmethod\tAID\n")
	else:
		open_out.write(utilities.PROTEIN_ID + "\tannotation\tfeature\tcategory\tmethod\tAID\n")
	outfile1 = re.sub(".tsv", ".attribute.tsv", outfile)
	open_out1 = open(outfile1, "w")
	open_out1.write("TID\t" + "AID" + "\t" + "key" + "\t" + "value" + "\n")
	mynum = 0
	
	for line in open_in:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^" + utilities.PROTEIN_FAMILY_ID, line) or re.search("^" + utilities.PROTEIN_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		if source == "protein_family":
			myclust = info[titles[utilities.PROTEIN_FAMILY_ID]]
		else:
			myclust = info[titles[utilities.PROTEIN_ID]]
		method = info[titles["method"]]
		study = info[titles["study"]]
		category = info[titles["category"]]
		mytype = info[titles["type"]]
		myid = myclust + "\t" + method + "\t" + mytype
		myclust_new = myclust 
		mynote = info[titles["note"]]
		myid1 = myclust + "\tDNA\tDNA_abundance"
		if not myid1 in annotation:
			if mynote == "good":
				mynote = "no_abundance"
			else:
				mynote = mynote + ";no_abundance"
		
		# cluster info
		if not myclust in flag_clust:
			flag_clust[myclust] = ""
			# note info: clarify the quality of the protein family: non-fungal Euk protein? unclassified_MSP (human-contaminated)?
			mystr = myclust_new + "\t" + mynote + "\t" + "quality" + "\t" + "note" + "\t" + "Quality_control\tNA"
			open_out.write(mystr + "\n")

			# cluster info
			mystr = myclust_new + "\t" + study + "\t" + "study" + "\t" + "project" + "\t" + "Shotgun\tNA"
			open_out.write(mystr + "\n")
			
			if not myclust in pep_cluster:
				# debug
				print("No cluster information!\t" + myclust)
				mystr = myclust_new + "\t" + "NA" + "\t" + "protein_family" + "\t" + "Denovo_clustering" + "\t" + "CD-hit\tNA"
				open_out.write(mystr + "\n")
			else:
				tmp1 = pep_cluster[myclust].split("\t")
				if tmp1[2] != "NA":
					myattr = myclust_new + "__" + "Denovo_clustering"
					mystr = myclust_new  + "\t"  + tmp1[0]+ "\t" + tmp1[1] + "\t" + "Denovo_clustering" + "\t" + "CD-hit\t" + myattr
					open_out.write(mystr + "\n")
					myline = myattr + "\n\n" + tmp1[2]
					mynum = output_attributes(myline, mynum, open_out1)
				else:
					mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + "Denovo_clustering" + "\t" + "CD-hit\t" + "NA"
					open_out.write(mystr + "\n")

			# mapping info
			if not myclust in mapping:
				# debug
				print("No taxonomy information!\t" + myclust)
				mystr = myclust_new + "\t" + "NA" + "\t" + "worse_homology" + "\t" + "UniRef90_homology" + "\t" + "UniRef90" + "\tNA"
				open_out.write(mystr + "\n")
			else:
				tmp1 = mapping[myclust].split("\t")
				if tmp1[2] != "NA":
					myattr = myclust_new + "__" + "UniRef90_homology"
					mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + "UniRef90_homology" + "\t" + "UniRef90" + "\t" + myattr
					open_out.write(mystr + "\n")
					myline = myattr + "\n\n" + tmp1[2]
					mynum = output_attributes(myline, mynum, open_out1)
				else:
					mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + "UniRef90_homology" + "\t" + "UniRef90" + "\t" + "NA"
					open_out.write(mystr + "\n")

			if not myclust in taxonomy:
				mystr = myclust_new + "\t" + "NA" + "\t" + "NA" + "\t" + "Taxonomy_characterization" + "\t" + "Taxonomy_annotation" + "\tNA"
				open_out.write(mystr + "\n")
			else:
				tmp1 = taxonomy[myclust].split("\t")
				if tmp1[2] != "NA":
					myattr = myclust_new + "__" + "Taxonomy_characterization"
					mystr = myclust_new + "\t"  + tmp1[0] + "\t" + tmp1[1] + "\t" + "Taxonomy_characterization" + "\t" + "Taxonomy_annotation" + "\t" + myattr
					open_out.write(mystr + "\n")
					myline = myattr + "\n\n" + tmp1[2]
					mynum = output_attributes(myline, mynum, open_out1)
				else:
					mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + "Taxonomy_characterization" + "\t" + "Taxonomy_annotation" + "\t" + "NA"
					open_out.write(mystr + "\n")

		# annotation info
		if myid in flag_type:
			continue
		flag_type[myid] = ""
		if myid in annotation:
			# uniref90 annotation
			if re.search("UniRef90", mytype):
				mysource = "UniRef90_characterization" + "\t" + "UniRef90"
				for tmp0 in annotation[myid]:
					tmp1 = tmp0.split("\t")
					if mytype == "UniRef90_unknown":
						mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + mysource + "\tNA"
						open_out.write(mystr + "\n")
					else:
						if tmp1[2] != "NA":
							myattr = myclust_new + "__" + tmp1[1]
							mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + mysource + "\t" + myattr
							open_out.write(mystr + "\n")
							myline = myattr + "\n\n" + tmp1[2]
							mynum = output_attributes(myline, mynum, open_out1)
						else:
							mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + mysource + "\t" + "NA"
							open_out.write(mystr + "\n")

			else: # Denovo annotation
				if mytype == "Denovo_signaling":
					method = "SignalP/Phobius"
				if mytype == "Denovo_transmembrane":
					method = "TMHMM/Phobius"
				mysource = "Denovo_characterization" + "\t" + method
				for tmp0 in annotation[myid]:
					tmp1 = tmp0.split("\t")
					if tmp1[2] != "NA":
						myattr = myclust_new + "__" + tmp1[1]
						mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + mysource + "\t" + myattr
						open_out.write(mystr + "\n")
						myline = myattr + "\n\n" + tmp1[2]
						mynum = output_attributes(myline, mynum, open_out1)
					else:
						mystr = myclust_new + "\t" + tmp1[0] + "\t" + tmp1[1] + "\t" + mysource + "\t" + "NA"
						open_out.write(mystr + "\n")
	# foreach line 
	open_out.close()

# finalize_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start finalize_annotation.py -l " + values.list + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get cluster info ......starting\n")
	pep_cluster = collect_cluster_info (config.protein_family, values.source)
	sys.stderr.write("Get cluster info ......done\n")
	
	### collect annotation info ###
	sys.stderr.write("Get annotation info ......starting\n")
	annotation = collect_annotation (values.list, values.source)
	taxonomy, mapping_tmp = collect_taxonomy_annotation (values.taxonomy)
	taxonomy_tmp, mapping = collect_taxonomy_annotation (values.uniref)
	sys.stderr.write("Get annotation info ......done\n")

	### finalize annotation for protein families ###
	sys.stderr.write("\nFinalize annotation for prioritization families ......starting\n")
	finalize_annotation (values.source, pep_cluster, annotation, taxonomy, mapping, values.annotation, values.output)
	sys.stderr.write("\nFinalize annotation for prioritization families ......done\n")

	sys.stderr.write("### Finish finalize_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

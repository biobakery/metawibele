#!/usr/bin/env python

"""
MetaWIBELE: mspminer_protein_family_taxonomy module
Summary taxonomy annotation for each protein family

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
from collections import Counter

try:
	from metawibele import config
	from metawibele import utilities
	from metawibele.common import utils
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Summary taxonomy annotation for each protein family
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--annotation",
	                    help='input the MSPs taxonomy annotation file',
	                    required=True)
	parser.add_argument('-s', "--source",
	                    help='specify the source of taxonomy for protein family',
	                    choices=["LCA", "Rep"],
	                    default="Rep")
	parser.add_argument('-o', "--output",
	                    help='output taxonomy annotation of protein family',
	                    required=True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_pep_cluster_info (clust_file):  # discovery_cohort.peptides.clust
	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	myclust_id = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			mym = re.search("cluster=([\d]+)", line)
			myclust_id = "Cluster_" + mym.group(1)
			if not myclust_id in cluster:
				cluster[myclust_id] = {}
			continue
		mym = re.search("([\S]+)", line)
		myid = mym.group(1)
		cluster[myclust_id][myid] = myclust
    # foreach line
	open_file.close()
	return cluster
# function collect_pep_cluster_info


#==============================================================
# collect taxonomy mapping info
#==============================================================
def collect_taxonomy_info (map_file):  
	taxa_map = {}
	titles = {}
	for line in utils.gzip_bzip2_biom_open_readlines (map_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^Taxon", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		# title line
		taxa_id = info[titles["Taxon"]]
		taxa_name = info[titles["Scientific_name"]]
		taxa_rank = info[titles["Rank"]]
		taxa_lineage = info[titles["Lineage"]]
		myid = taxa_name
		tmp = taxa_lineage.split("|")
		if re.search("__", tmp[-1]):
			mym = re.search("([^\_]+)__([\S]+)", tmp[-1])
			myrank = mym.group(1)
			myid = mym.group(2)
			if myrank == "k":
				myrank = "Kingdom"
			if myrank == "p":
				myrank = "Phylum"
			if myrank == "c":
				myrank = "Class"
			if myrank == "o":
				myrank = "Order"
			if myrank == "f":
				myrank = "Family"
			if myrank == "g":
				myrank = "Genus"
			if myrank == "s":
				myrank = "Species"
			if myrank == "t":
				myrank = "Terminal"
			taxa_rank = myrank
		taxa_map[myid] = taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage
	# foreach line
	
	return taxa_map
# collect_taxonomy_info


#==============================================================
# extract taxon of MSP-based annotation for each gene
#==============================================================
def extract_taxon_info (msp_file):  
	# collect info
	genes = {}
	taxa = {}
	titles = {}
	ann_title = ""
	open_file = open(msp_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search(utilities.PROTEIN_ID, line):
			for item in info:
				titles[item] = info.index(item)
			ann_title = "\t".join(info[1:titles["taxa_id"]]) + "\tnote\tmsp_name\tmsp_taxa_name\tmsp_taxa_id"
			continue
		myid = info[titles[utilities.PROTEIN_ID]]
		taxID = info[titles["taxa_id"]]
		taxa_name = info[titles["taxa_name"]]
		taxa_rank = info[titles["taxa_rank"]]
		tmp = info[titles["taxa_lineage"]].split("|")
		taxa_lineage = "NA"
		for item in tmp:
			if re.search("__", item):
				if taxa_lineage == "NA":
					taxa_lineage = item
				else:
					taxa_lineage = taxa_lineage + "|" + item
		#organism = info[titles["organism"]]
		map_type = info[titles["map_type"]]
		genes[myid] = "\t".join(info[1:titles["taxa_id"]]) + "\t" + info[titles["note"]] + "\t" + info[titles["msp_name"]] + "\t" + info[titles["msp_taxa_name"]] + "\t" + info[titles["msp_taxa_id"]] + "\n" + taxa_lineage
		if taxa_lineage != "NA":
			tmp = taxa_lineage.split("|")
			for item in tmp:
				if re.search("__", item):
					mym = re.search("^([^\_]+)__([\S]+)", item)
					myrank = mym.group(1)
					myname = mym.group(2)
					if myname == "NA":
						#continue
						myrank = "Unclassified"
					if myrank == "k":
						myrank = "Kingdom"
					if myrank == "p":
						myrank = "Phylum"
					if myrank == "c":
						myrank = "Class"
					if myrank == "o":
						myrank = "Order"
					if myrank == "f":
						myrank = "Family"
					if myrank == "g":
						myrank = "Genus"
					if myrank == "s":
						myrank = "Species"
					if myrank == "t":
						myrank = "Terminal"
					#myname = re.sub("_", " ", myname)
					if not myrank in taxa:
						taxa[myrank] = {}
					#if not gene_name in taxa[myrank]:
					taxa[myrank][myid] = myname
				# if taxon
			# foreach taxon rank
		# if taxon exists
	# foreach gene
	open_file.close()
	return genes, taxa, ann_title
# extract_taxon_info


#==============================================================
# taxonomy annotation for protein family
#==============================================================
def taxonomy_annotation (genes, taxa, ann_title, min_cutoff, taxa_map, taxa_source, pep_cluster, outfile):
	outfile2 = re.sub(".tsv", ".all.tsv", outfile)
	outfile3 = re.sub("_proteinfamilies_", "_protein_", outfile)
	outfile4 = re.sub("_proteinfamilies_", "_protein_", outfile2)
	taxa_level = ["Terminal", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Unclassified"]
	gene_taxa = {}
	cluster_taxa = {}
	rep_taxa = {}
	rep_msp = {}
	for myclust in sorted(pep_cluster.keys()):
		taxa_flag = 0
		taxa_name = "NA"
		taxa_num = {}
		total = 0
		for myid in sorted(pep_cluster[myclust].keys()):	
			myrep = pep_cluster[myclust][myid]
			if myid in genes:
				myinfo1, myline = genes[myid].split("\n")
				myinfo2 = "NA\tNA\tNA\tNA\tNA\tNA\tUnclassified\tNA"
				tmp = myline.split("|")
				if re.search("__", tmp[-1]):
					mym = re.search("__([\S]+)", tmp[-1])
					mytaxa1 = mym.group(1)
					if mytaxa1 in taxa_map:
						taxa_id, taxa_name, taxa_rank, taxa_lineage = taxa_map[mytaxa1].split("\t")
						myinfo2 = taxa_id + "\t" + taxa_name + "\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage
				gene_taxa[myid] = myinfo1 + "\t" + myinfo2
				if myid == myrep:
					tmp2 = myinfo2.split("\t")
					rep_taxa[myclust] = myinfo1 + "\n" + tmp2[0] + "\t" + tmp2[1]
					tmp3 = myinfo1.split("\t")
					rep_msp[myclust] = tmp2[-2] + "\t" + tmp3[-3]
					# debug
					#print(myclust + "\t" + tmp2[-2] + "\t" + tmp3[-1])
				total = total + 1
				for mytaxa in taxa.keys():
					if myid in taxa[mytaxa]:
						myname = taxa[mytaxa][myid]
						if not mytaxa in taxa_num:
							taxa_num[mytaxa] = []
						taxa_num[mytaxa].append(myname)
				# foreach taxon rank
			# if myid
		# foreach gene in the cluster
		if total == 0:
			taxa_flag = 1
			taxa_name = "NA"
		for mylevel in taxa_level:
			if taxa_flag == 1:
				break
			if mylevel in taxa_num:
				tmp = Counter(taxa_num[mylevel])
				first_tax = "NA"
				first_num = 0
				try:
					for (mytax, tax_num) in sorted(tmp.iteritems(), key = lambda d:d[1], reverse = True): # python2
						first_tax = mytax
						first_num = int(tax_num)
						break	
				except:
					for (mytax, tax_num) in sorted(tmp.items(), key = lambda d:d[1], reverse = True): # python3
						first_tax = mytax
						first_num = int(tax_num)
						break	
				# foreach type of taxon
				mydiff = first_num
				myper = float(mydiff) / float(total)
				if myper >= float(min_cutoff):	# above the consistency cutoff 
					taxa_flag = 1
					taxa_name = first_tax
			# if taxa_num
		# foreach taxon level
		if taxa_flag == 0:	# still unclassified
			taxa_name = "NA"
		if myclust in rep_taxa:
			myinfo1, tmp = rep_taxa[myclust].split("\n")
			myrep_id, myrep = tmp.split("\t")
			mylac = taxa_name
			mylac_id = "NA"
			if taxa_name in taxa_map:
				tmp = taxa_map[taxa_name].split("\t")
				mylac_id = tmp[0]
				mylac = tmp[1]
			myname = taxa_name
			if taxa_source == "Rep":
				myname = myrep
				myname = re.sub(" ", "_", myname)
			taxa_info = "NA\tNA\tUnclassified\tNA"
			if myname in taxa_map:
				taxa_info = taxa_map[myname]
				tmp = taxa_info.split("\t")
				if re.search("Superkingdom", tmp[2]) or re.search("Kingdom", tmp[2]): # kindom as unclassified
					taxa_info = "NA\tNA\tUnclassified\tNA"
			if not myclust in cluster_taxa:
				cluster_taxa[myclust] = myinfo1 + "\t" + mylac + "\t" + mylac_id + "\t" + myrep + "\t" + myrep_id + "\t" + taxa_info
		# if rep_taxa 
	# foreach cluster

	# output info
	open_out1 = open(outfile, "w")
	open_out2 = open(outfile2, "w")
	open_out3 = open(outfile3, "w")
	open_out4 = open(outfile4, "w")
	open_out1.write(utilities.PROTEIN_FAMILY_ID + "\t"  + ann_title + "\t" + "MSP_Tax\tMSP_TaxID\tMSP_Rep_Tax\tMSP_Rep_TaxID\t" + "taxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\n")
	open_out2.write(utilities.PROTEIN_FAMILY_ID + "\t"  + ann_title + "\t" + "MSP_Tax\tMSP_TaxID\tMSP_Rep_Tax\tMSP_Rep_TaxID\t" + "taxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\n")
	open_out3.write(utilities.PROTEIN_ID + "\t"  + ann_title + "\t" + "MSP_Tax\tMSP_TaxID\tMSP_Rep_Tax\tMSP_Rep_TaxID\t" + "taxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\n")
	open_out4.write(utilities.PROTEIN_ID + "\t"  + ann_title + "\t" + "MSP_Tax\tMSP_TaxID\tMSP_Rep_Tax\tMSP_Rep_TaxID\t" + "taxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\n")
	taxa_level = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Terminal"]
	for myclust in sorted(cluster_taxa.keys()):
		info = cluster_taxa[myclust].split("\t")
		mymsp = info[-11]
		myrank = info[-2]
		if myclust in rep_msp:
			myrank, mymsp = rep_msp[myclust].split("\t")
		if mymsp != "msp_unknown" and myrank == "Unclassified":
			mynote = "unclassified_MSP"
			info[-12] = info[-12] + ";" + mynote
		tmp_str = "\t".join(info)
		open_out1.write(myclust + "\t" + tmp_str + "\n")
		#info = cluster_taxa[myclust].split("\t")
		pre_line = myclust + "\t" + "\t".join(info[0:len(info)-4])
		myinfo = info[-1].split("|")
		mylevel = {}
		for item in myinfo:
			mytaxa= "NA\tNA\tUnclassified\tNA"
			if re.search("__", item):
				mym = re.search("__([\S]+)", item)
				myname = mym.group(1)
				if myname in taxa_map:
					mytaxa = taxa_map[myname]
					mytmp = mytaxa.split("\t")
					taxa_id = mytmp[0]
					taxa_name = mytmp[1]
					taxa_rank = mytmp[2]
					open_out2.write(pre_line + "\t" + mytaxa + "\n")
					mylevel[taxa_rank] = ""
			#if item == "NA":
			#	open_out2.write(pre_line + "\t" + mytaxa + "\n") 
		# foreach item
		for myl in taxa_level:
			if myl in mylevel:
				continue
			else:
				taxa_name = "Unclassified"
				taxa_rank = myl
				taxa_lineage = "NA"
				taxa_id = "NA"
				open_out2.write(pre_line + "\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage + "\n")
		# foreach taxonomy level
	# foreach cluster
	open_out1.close()
	open_out2.close()
	for myid in sorted(gene_taxa.keys()):
		info = gene_taxa[myid].split("\t")
		mymsp = info[-11]
		myrank = info[-2]
		if mymsp != "msp_unknown" and myrank == "Unclassified":
			mynote = "unclassified_MSP"
			info[-12] = info[-12] + ";" + mynote
		tmp_str = "\t".join(info)
		open_out3.write(myid + "\t" + tmp_str + "\n")
		#info = gene_taxa[myid].split("\t")
		pre_line = myid + "\t" + "\t".join(info[0:len(info)-4])
		myinfo = info[-1].split("|")
		mylevel = {}
		for item in myinfo:
			mytaxa= "NA\tNA\tUnclassified\tNA"
			if re.search("__", item):
				mym = re.search("__([\S]+)", item)
				myname = mym.group(1)
				if myname in taxa_map:
					mytaxa = taxa_map[myname]
					open_out4.write(pre_line + "\t" + mytaxa + "\n")
					mytmp = mytaxa.split("\t")
					taxa_rank = mytmp[2]
					mylevel[taxa_rank] = ""
			#if item == "NA":
			#	open_out4.write(pre_line + "\t" + mytaxa + "\n") 
		# foreach item
		for myl in taxa_level:
			if myl in mylevel:
				continue
			else:
				taxa_name = "Unclassified"
				taxa_rank = myl
				taxa_lineage = "NA"
				taxa_id = "NA"
				open_out4.write(pre_line + "\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage + "\n")
		# foreach taxonomy level
	# foreach cluster
	open_out3.close()
	open_out4.close()
# taxonomy_annotation 


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	sys.stderr.write("### Start mspminer_protein_family_taxonomy.py -a " + values.annotation + " ####\n")
	

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	pep_cluster = collect_pep_cluster_info (config.protein_family)
	taxa_map = collect_taxonomy_info (config.taxonomy_database)
	genes, taxa, ann_title = extract_taxon_info (values.annotation)
	taxonomy_annotation (genes, taxa, ann_title, config.tshld_lca, taxa_map, values.source, pep_cluster, values.output)
	sys.stderr.write("Get info ......done\n")

	sys.stderr.write("### Finish mspminer_protein_family_taxonomy.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

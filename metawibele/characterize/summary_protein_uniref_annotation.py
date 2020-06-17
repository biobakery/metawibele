#!/usr/bin/env python

"""
MetaWIBELE: summary_protein_uniref_annotation module
Summary UniRef annotation for each ORF

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
	from metawibele.common import utils
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Summary UniRef annotation for each ORF
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-a', "--annotation",
	                    help='input uniref functional annotation info file',
	                    required=True)
	parser.add_argument('-m', "--mapping",
	                    help='input UniRef mapping summary file',
	                    required=True)
	parser.add_argument('-t', "--type",
	                    help='specify the taxonomy source of UniRef',
						choices = ["LCA", "Rep"],
	                    required=True,
	                    default="Rep")
	parser.add_argument('-o', "--output",
	                    help='output cluster distribution file',
	                    required=True)
	values = parser.parse_args()

	return values
# get_args


#==============================================================
# collect taxonomy info for taxaID
#==============================================================
def collect_taxonomy_info (taxa_file, taxa_hits):  # uniprot_taxonomy.map.tsv
	taxa = {}
	taxa_map = {}
	titles = {}
	for line in utils.gzip_bzip2_biom_open_readlines (taxa_file): 
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^Taxon", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		# title line
		mytaxa = info[titles["Taxon"]]
		if not mytaxa in taxa_hits:
			continue
		myname = info[titles["Scientific_name"]]
		myline = info[titles["Lineage"]]
		myrank = info[titles["Rank"]]
		myinfo = myline.split("|")
		if re.search("unclassified_sequences", myinfo[0]):
			continue
		myid = myname
		if re.search("__", myinfo[-1]):
			mym = re.search("([^\_]+)__([\S]+)", myinfo[-1])
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
		taxa_map[myid] = mytaxa + "\t" + myname + "\t" + myrank + "\t" + myline
		taxa[mytaxa] = myname + "\t" + myrank + "\t" + myline
	# foreach line
	
	return taxa, taxa_map
# function collect_taxonomy_info


#==============================================================
# collect taxonomy info for UniRefID
#==============================================================
def collect_uniref_taxonomy_info (uniref_file, spe_type, hits):  # uniref90.ann.all.tsv 
	uniref_taxa = {}
	taxa_hits = {}
	titles = {}
	open_file = open(uniref_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^ID", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		uniref_id = info[titles["ID"]]
		if not uniref_id in hits:
			continue
		taxa_id = info[titles["TaxID"]]
		taxa = info[titles["Tax"]]
		reptaxa_id = info[titles["Rep_TaxID"]]
		reptaxa = info[titles["Rep_Tax"]]
		if re.search("_UPI", uniref_id) and reptaxa_id == "NA":	# use LCA as representatives for UniParc item, e.g. UniRef90_UPI000682ABCC
			reptaxa_id = taxa_id
			reptaxa = taxa
		mytaxa_id = reptaxa_id
		if spe_type == "LCA":
			mytaxa_id = taxa_id
		if spe_type == "Rep":
			mytaxa_id = reptaxa_id
		taxa_hits[mytaxa_id] = ""
		detail = info[titles["Protein_names"]]
		#org = info[titles["Organism"]]
		uniprot = info[titles["UniProtKB"]]
		uniref_taxa[uniref_id] = mytaxa_id + "\n" + detail + "\t" + taxa + "\t" + taxa_id + "\t" + reptaxa + "\t" + reptaxa_id + "\t" + uniprot
	# foreach line
	open_file.close()

	return uniref_taxa, taxa_hits
# collect_uniref_taxonomy_info


#==============================================================
# collect cluster info
#==============================================================
def collect_pep_cluster_info (clust_file):	# discovery_cohort.peptides.clust
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
			#if not myclust_id in cluster:
			#	cluster[myclust_id] = myclust
			continue
		mym = re.search("([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myid
	# foreach line
	open_file.close()
	return cluster
# function collect_pep_cluster_info

def collect_gene_cluster_info (clust_file):	# discovery_cohort.genes.clust
	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myclust
	# foreach line
	open_file.close()
	return cluster
# function collect_gene_cluster_info


#==============================================================
# collect UniRef summary info
#==============================================================
def collect_uniref_mapping_info (uniref_stat, uniref_taxa, taxa, taxa_map):	
	uniref = {}
	titles = {}
	taxa_level = ['g', 'f', 'o', 'c', 'p']
	open_file = open(uniref_stat, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search(utilities.PROTEIN_ID, line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myid = info[titles[utilities.PROTEIN_ID]]
		map_type1 = info[titles["query_type"]]
		map_type2 = info[titles["mutual_type"]]
		uniref_id = info[titles["subject"]]
		map_iden = info[titles["identity"]]
		map_cov1 = info[titles["query_coverage"]]
		map_cov2 = info[titles["mutual_coverage"]]
		mytaxa = "NA\nNA\tNA\tNA\tNA\tNA\tNA\tNA"
		if uniref_id in uniref_taxa:
			mytaxa = uniref_taxa[uniref_id]

		"""
			if map_type2 == "low_confidence" or map_type1 == "low_confidence": # weak homology
				myflag = "NA"
				if float(map_iden) >= 85 and float(map_cov1) >= 0.8: # genus level
					myflag = "g"
				else:
					if float(map_iden) < 85 and float(map_iden) >= 65 and float(map_cov1) >= 0.8: # phylum level
						myflag = "p"
				mytax_id, tmp = mytaxa.split("\n")
				if myflag == "NA":
					mytaxa = "NA" + "\n" + tmp
				else:
					myflag2 = 0
					if mytax_id in taxa:
						tmp1 = taxa[mytax_id].split("\t")
						tmp2 = tmp1[-1].split("|")
						tmp2.reverse()
						myindex = taxa_level.index(myflag)
						myflag1 = 0
						while myindex < len(taxa_level):
							myitem = taxa_level[myindex]
							if re.search(myitem + "__", tmp1[-1]):  # find the level
								myflag = myitem
								myflag1= 1
								break
							myindex = myindex + 1
						if myflag1 == 1:
							for item in tmp2:
								if re.search("__", item):
									mym = re.search("([^\_]+)__([\S]+)", item)
									tmp_myrank = mym.group(1)
									tmp_myid = mym.group(2)
									if tmp_myrank == myflag:
										if tmp_myid in taxa_map:
											tmp3 = taxa_map[tmp_myid].split("\t")
											new_id = tmp3[0]
											mytaxa = new_id + "\n" + tmp
											myflag2 = 1
											break
					if myflag2 == 0:
						mytaxa = "NA" + "\n" + tmp
		"""

		uniref[myid] = map_type1 + "\t" + map_type2 + "\t" + map_iden + "\t" + map_cov1 + "\t" + map_cov2 + "\n\n" + uniref_id + "\n\n" + mytaxa
	# foreach line
	open_file.close()
	return uniref
# collect_uniref_mapping_info


#==============================================================
# collect pre-annotation info
#==============================================================
def collect_ann_info (ann_file, spe_type):	# summary_peptide_family_annotation.tsv  
	anns = {}
	ann_type = {}
	hits = {}
	uniref_taxa = {}
	taxa_hits = {}
	titles = {}
	open_file = open(ann_file, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	myindex = 0
	while myindex < len(info):
		titles[info[myindex]] = myindex
		myindex = myindex + 1
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^" + utilities.PROTEIN_ID, line):
			continue
		info = line.split("\t")
		myid = info[titles[utilities.PROTEIN_ID]]
		mytype = info[titles["type"]]
		detail = info[titles["Protein_names"]]
		mytax = info[titles["Tax"]]
		mytaxID = info[titles["TaxID"]]
		myreptax = info[titles["Rep_Tax"]]
		myreptaxID = info[titles["Rep_TaxID"]]
		#myorg = info[titles["organism"]]
		myuniprot = info[titles["UniProtKB"]]
		myunirefID = info[titles["unirefID"]]
		hits[myunirefID] = ""
		mystr = detail + "\t" + mytax + "\t" + mytaxID + "\t" + myreptax + "\t" + myreptaxID + "\t" + myuniprot
		if not re.search("^UniRef", mytype):
			continue
		if not myid in ann_type:
			ann_type[myid] = {}
		ann_type[myid][mytype] = ""
		anns[myid] = mystr
		
		if re.search("_UPI", myunirefID) and myreptaxID == "NA":	# use LCA as representatives for UniParc item, e.g. UniRef90_UPI000682ABCC
			myreptaxID = mytaxID
			myreptax = mytax
		mytaxa_id = myreptaxID
		if spe_type == "LCA" or spe_type == "lca":
			mytaxa_id = mytaxID
		taxa_hits[mytaxa_id] = ""
		uniref_taxa[myunirefID] = mytaxa_id + "\n" + detail + "\t" + mytax + "\t" + mytaxID + "\t" + myreptax + "\t" + myreptaxID + "\t" + myuniprot
	# foreach line
	open_file.close()
	
	return anns, ann_type, hits, uniref_taxa, taxa_hits
# collect_anns_info


#==============================================================
# assign uniref flag to peptide families 
#==============================================================
def assign_annotation (identity_cutoff, coverage_cutoff, cluster, uniref, anns, ann_type, taxa, taxa_map, study, outfile):
	outfile2 = re.sub(".tsv", ".all.tsv", outfile)
	open_out = open(outfile, "w")
	open_out2 = open(outfile2, "w")
	open_out.write(utilities.PROTEIN_ID + "\tstudy\tmap_type\tquery_type\tmutual_type\tidentity\tquery_coverage\tmutual_coverage\tdetail\tTax\tTaxID\tRep_Tax\tRep_TaxID\tUniProtKB\tunirefID\ttaxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\tnote\n")
	open_out2.write(utilities.PROTEIN_ID + "\tstudy\tmap_type\tquery_type\tmutual_type\tidentity\tquery_coverage\tmutual_coverage\tdetail\tTax\tTaxID\tRep_Tax\tRep_TaxID\tUniProtKB\tunirefID\ttaxa_id\ttaxa_name\ttaxa_rank\ttaxa_lineage\tnote\n")
	taxa_level = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Terminal"]
	for myclust in sorted(anns.keys()): # foreach peptide family
		mynote = "good"
		if myclust in cluster:
			myid = cluster[myclust]	 # protein id
			#if not myid in gene_cluster:  # no corresponding gene cluster
			#	print("Peptide ID has no corresponding gene cluster!\t" + myid)
			#	continue
			#gene_id = gene_cluster[myid]
			gene_id = myid
			if gene_id in uniref:
				mytype, uniref_id, uniref_taxa = uniref[gene_id].split("\n\n")
				myquery, mymutual, myiden, query_cov, mutual_cov = mytype.split("\t")
				mytaxa, mystr = uniref_taxa.split("\n")
				taxa_line = "NA\tUnclassified\tNA"
				if mytaxa in taxa:
					taxa_line = taxa[mytaxa]
				if re.search("k__Eukaryota", taxa_line) and not re.search("k__Fungi", taxa_line):
						# debug
						#print("Filter out Eukaryota proteins\t" + myclust + "\t" + uniref_id + "\t" + taxa_line)
						#continue
					#if float(myiden) >= 50 and float(mutual_cov) >= 0.8:
					mynote = "non-fungal Eukaryota proteins"
				if mymutual == "high_confidence":
					mytype = "UniRef90_characterized"
				if mymutual == "low_confidence":
					mytype = "UniRef90_weak_homology"
					if float(myiden) < float(identity_cutoff) or float(mutual_cov) < float(coverage_cutoff):
						mytype = "UniRef90_worse_homology"
				if mymutual == "no_hit":
					mytype = "UniRef90_worse_homology"
				# correct
				if not "UniRef90_unknown" in ann_type[myclust]:
					mytype = "UniRef90_characterized"
					if "UniRef90_uncharacterized" in ann_type[myclust]:
						# debug
						#print("Correct uncharacterized type!\t" + myclust)
						mytype = "UniRef90_uncharacterized"
				if "UniRef90_unknown" in ann_type[myclust]:
					if mytype == "UniRef90_characterized":
						# debug
						print("Change characterized type to unknown!\t" + myclust)
						mytype = "UniRef90_weak_homology"
					if mytype == "UniRef90_uncharacterized":
						# debug
						print("Change uncharacterized type to unknown!\t" + myclust)
						mytype = "UniRef90_weak_homology"
				if mytype == "UniRef90_weak_homology":
					anns[myclust] = mystr
				if mytype == "UniRef90_worse_homology":
					uniref_id = "NA"
					taxa_line = "NA\tUnclassified\tNA"
					mynote = "good"
				open_out.write(myclust + "\t" + study + "\t" + mytype + "\t" + myquery + "\t" + mymutual + "\t" + myiden + "\t" + query_cov + "\t" + mutual_cov + "\t" + anns[myclust] + "\t" + uniref_id + "\t" + mytaxa + "\t" + taxa_line + "\t" + mynote + "\n")
				# output all taxonomy level
				pre_line = myclust + "\t" + study + "\t" + mytype + "\t" + myquery + "\t" + mymutual + "\t" + myiden + "\t" + query_cov + "\t" + mutual_cov + "\t" + anns[myclust] + "\t" + uniref_id
				tmp = taxa_line.split("\t")
				myinfo = tmp[-1].split("|")
				mylevel = {}
				#item_sub = "\t".join(tmp[0:len(tmp)-3])
				for item in myinfo:
					myrank = "NA"
					myname = "NA"
					mytmp = "NA\tNA\tNA\tNA"
					if re.search("k__", item):
						mym = re.search("k__([\S]+)", item)
						myrank = "Kingdom"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if re.search("p__", item):
						mym = re.search("p__([\S]+)", item)
						myrank = "Phylum"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if re.search("c__", item):
						mym = re.search("c__([\S]+)", item)
						myrank = "Class"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if re.search("o__", item):
						mym = re.search("o__([\S]+)", item)
						myrank = "Order"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if re.search("f__", item):
						mym = re.search("f__([\S]+)", item)
						myrank = "Family"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if re.search("g__", item):
						mym = re.search("g__([\S]+)", item)
						myrank = "Genus"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if re.search("s__", item):
						mym = re.search("s__([\S]+)", item)
						myrank = "Species"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if re.search("t__", item):
						mym = re.search("t__([\S]+)", item)
						myrank = "Terminal"
						#myname = re.sub("_", " ", mym.group(1))
						myname = mym.group(1)
						if myname in taxa_map:
							mytmp = taxa_map[myname]
						mylevel[myrank] = ""
					if myrank != "NA" and myname != "NA":
						#myname = re.sub("_", " ", myname)
						open_out2.write(pre_line + "\t" + mytmp + "\t"+ mynote + "\n")
					#if item == "NA":
					#	myname = tmp[-3]
					#	myrank = tmp[-2]
					#	open_out2.write(pre_line + "\t" + myname + "\t" + myrank + "\t" + tmp[-1] + "\n")
				# foreach level
				for myl in taxa_level:
					if myl in mylevel:
						continue
					else:
						myname = "Unclassified"
						myrank = myl
						taxa_lineage = "NA"
						taxa_id = "NA"
						taxa_name = "Unclassified"
						taxa_rank = myl
						open_out2.write(pre_line + "\t" + taxa_id + "\t" + taxa_name + "\t" + taxa_rank + "\t" + taxa_lineage + "\t" + mynote + "\n")
			else:
				print("No uniref ID info!\t" + myclust)
		# cluster name
		else:
			print("No cluster info!\t" + myclust)
	# foreach peptide cluster
	open_out.close()	
	open_out2.close()	
# assign_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()


	sys.stderr.write("### Start summary_protein_uniref_annotation.py -a " + values.annotation + " ####\n")

	### collect cluster info ###
	sys.stderr.write("Get info ......starting\n")
	pep_cluster = collect_pep_cluster_info (config.protein_family)
	anns, ann_type, hits, uniref_taxa, taxa_hits = collect_ann_info (values.annotation, values.type)
	taxa, taxa_map = collect_taxonomy_info (config.taxonomy_database, taxa_hits)
	uniref = collect_uniref_mapping_info (values.mapping, uniref_taxa, taxa, taxa_map)
	sys.stderr.write("Get info ......done\n")

	### assign annotation to peptide families ###
	sys.stderr.write("\nAssign annotation to peptide families ......starting\n")
	assign_annotation (config.tshld_identity, config.tshld_coverage, pep_cluster, uniref, anns, ann_type, taxa, taxa_map, config.study, values.output)
	sys.stderr.write("\nAssign annotation to peptide families ......done\n")

	sys.stderr.write("### Finish summary_protein_uniref_annotation.py ####\n\n\n")

# end: main

if __name__ == '__main__':
	main()

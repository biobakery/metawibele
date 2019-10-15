#!/usr/bin/env python

"""
MeteWIBELE workflow: MeteWIBELE characterization workflow
1) clustering into protein families
2) homology-based annotation
3) sequence-based annotation
4) abundance-based annotation
5) combination

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
import os, fnmatch
import logging

# import the workflow class from anadama2
from anadama2 import Workflow

# import the library of MetaWIBELE tasks for characterization
from metawibele.tasks import characterization

# import the utilities functions and config settings from MetaWIBELE
from metawibele import utilities, config


def parse_cli_arguments ():
	'''
	 Parses any command-line arguments passed into the workflow.
	the version number will appear when running this script with the "--version" option
	the description will appear when running this script with the "--help" option
	create a workflow instance, providing the version number and description
	'''

	workflow = Workflow(version = "0.0.1", description = "A workflow for MetaWIBELE characterization")

	# add the custom arguments to the workflow
	workflow.add_argument("threads",
	                      desc = "number of threads/cores for each task to use",
	                      default = 20)
	workflow.add_argument("config",
	                      desc = "the configuration file of characterization analysis",
	                      default = "characterization.cfg")
	workflow.add_argument("clustering",
	                      desc = "indicates whether or not cluster proteins into protein families",
	                      default = True,
	                      action = "store_true")
	workflow.add_argument("homology-based",
	                      desc = "indicates whether or not annotate protein families based on homology information",
	                      default = True,
	                      action = "store_true")
	workflow.add_argument("sequence-based",
	                      desc = "indicates whether or not annotate protein families based on sequence information",
	                      default = True,
	                      action = "store_true")
	workflow.add_argument("abundance-based",
	                      desc = "indicates whether or not annotate protein families based on abundance information",
	                      default = True,
	                      action = "store_true")
	workflow.add_argument("split-number",
	                      desc="indicates number of spliting files for annotation based on sequence information",
	                      default = 20)
	workflow.add_argument("rna-count",
	                      desc = "indicates RNA count tables if is available",
	                      default = "none")
	workflow.add_argument("finalized",
	                      desc = "indicates whether or not finalize annotations for protein families",
	                      default = True,
	                      action = "store_true")

	return workflow



def get_method_config (config_file):
	'''
	:param config_file:
	:return: config info for each analysis
	'''
	print('get_method_config')

	config_items = config.read_user_edit_config_file(config_file)
	homology_conf = {}
	sequence_conf = {}
	abundance_conf = {}
	final_conf = {}
	values = ["yes", "no", "Yes", "No"]

	if "homology_based" in config_items:
		for name in config_items["homology_based"].keys():
			myvalue = config_items["homology_based"][name]
			if not myvalue in values:
				print('The config value can not be recognized. Please check your config file!')
				continue
			homology_conf[name] = myvalue
		# for each method

	if "sequence_based" in config_items:
		for name in config_items["sequence_based"].keys():
			myvalue = config_items["sequence_based"][name]
			if not myvalue in values:
				print('The config value can not be recognized. Please check your config file!')
				continue
			sequence_conf[name] = myvalue
		# for each method

	if "abundance_based" in config_items:
		for name in config_items["abundance_based"].keys():
			myvalue = config_items["abundance_based"][name]
			abundance_conf[name] = myvalue
		# for each method

	if "finalized" in config_items:
		for name in config_items["finalized"].keys():
			myvalue = config_items["finalized"][name]
			if not myvalue in values:
				print('The config value can not be recognized. Please check your config file!')
				continue
			final_conf[name] = myvalue
		# for each method

	return homology_conf, sequence_conf, abundance_conf, final_conf



def main(workflow):
	'''
	Characterization of protein families
	'''

	# get arguments
	args = workflow.parse_args()

	# get configuration info
	homology_conf, sequence_conf, abundance_conf, final_conf = get_method_config(args.config)

	# input and output folder
	input_dir = args.input
	output_dir = config.annotation_dir

	# get all input files
	gene_catalog = config.gene_catalog
	gene_catalog_seq = config.gene_catalog_prot
	gene_catalog_count = config.gene_catalog_count
	gene_catalog_rna_count = args.rna_count

	# get all output files
	protein_family = config.protein_family
	protein_family_seq = config.protein_family_prot_seq
	protein_family_relab = config.protein_family_relab
	protein_family_ann = config.protein_family_ann
	protein_family_attr = config.protein_family_attr
	protein_family_ann_list = {}
	protein_ann_list = {}
	uniref_taxonomy_family = output_dir + "/" + config.basename + "_proteinfamilies_annotation.uniref90_annotation.tsv"
	uniref_taxonomy = output_dir + "/" + config.basename + "_protein_annotation.uniref90_annotation.tsv"
	taxonomy_annotation_family = output_dir + "/" + config.basename + "_proteinfamilies_annotation.MSPminer_taxonomy.tsv"
	taxonomy_annotation = output_dir + "/" + config.basename + "_protein_annotation.MSPminer_taxonomy.tsv"


	### STEP #1: clustering ###
	# if clustering action is provided, then cluster proteins into protein families
	if args.clustering:
		myprotein_family, myprotein_family_output_folder = characterization.clustering (workflow, gene_catalog_seq,
		                                                                                args.threads,
		                                                                                output_dir, protein_family, protein_family_seq)

	### STEP #2: homology-based annotation ###
	# if homology-based action is provided, then directly extract annotations from database (UniProt)
	if args.homology_based:
		myprotein_family_ann, myprotein_ann, homology_output_folder = characterization.homology_based_annotation (workflow, homology_conf,
		                                                                                                          gene_catalog_seq,
		                                                                                                          args.threads,
		                                                                                                          output_dir, uniref_taxonomy_family, uniref_taxonomy,
		                                                                                                          protein_family_ann_list, protein_ann_list)

	### STEP #3: sequence-based annotation ###
	# if sequence-based action is provided, then do annotations based on sequence information
	if args.sequence_based:
		myprotein_family_ann, myprotein_ann, sequence_output_folder = characterization.sequence_based_annotation (workflow, sequence_conf,
		                                                                                                          gene_catalog_seq,
		                                                                                                          args.split_number, args.threads,
		                                                                                                          output_dir, protein_family_ann_list, protein_ann_list)


	### STEP #4: abundance-based annotation ###
	# if abundance-based action is provided, then do annotations based on abundance information
	if args.abundance_based:
		myprotein_family_ann, myprotein_ann, abundance_output_folder = characterization.abundance_based_annotation (workflow, abundance_conf,
		                                                                                                            gene_catalog, gene_catalog_count, gene_catalog_rna_count,
		                                                                                                            uniref_annotation_family, uniref_taxonomy,
		                                                                                                            args.split_number, args.threads,
		                                                                                                            output_dir, protein_family_relab, taxonomy_annotation_family, taxonomy_annotation,
		                                                                                                            protein_family_ann_list, protein_ann_list)

	### STEP #4: finalize annotations ###
	if args.finalized:
		myprotein_family_ann, myprotein_family_attr, annotation_output_folder = characterization.finalized_annotation (workflow, final_conf,
		                                                                                                               protein_family_ann_list, protein_ann_list,
		                                                                                                               protein_family_ann_list_file, protein_ann_list_file,
		                                                                                                               uniref_annotation_family, uniref_annotation,
		                                                                                                               taxonomy_annotation_family, taxonomy_annotation,
		                                                                                                               args.threads,
		                                                                                                               output_dir, protein_family_ann, protein_family_attr)

	### start the workflow
	workflow.go()



if __name__ == "__main__":
	main(parse_cli_arguments())



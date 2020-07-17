#!/usr/bin/env python

"""
MeteWIBELE workflow: MeteWIBELE characterization workflow
1) clustering into protein families
2) protein family annotation
3) domain/motif annotation
4) abundance annotation
5) integration

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
try:
	from metawibele.tasks import characterization
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
		         " Please check your install.")

# import the utilities functions and config settings from MetaWIBELE
try:
	from metawibele import utilities, config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
		         " Please check your install.")


VERSION = config.version

def parse_cli_arguments ():
	'''
	 Parses any command-line arguments passed into the workflow.
	the version number will appear when running this script with the "--version" option
	the description will appear when running this script with the "--help" option
	create a workflow instance, providing the version number and description
	'''

	workflow = Workflow(version = VERSION, description = "A workflow for MetaWIBELE characterization", remove_options=["input","output"])

	# add the custom arguments to the workflow
	workflow.add_argument("threads",
	                      desc = "number of threads/cores for each task to use",
	                      default = None)
	workflow.add_argument("characterization-config",
	                      desc = "the configuration file of characterization analysis",
	                      default = None)
	workflow.add_argument("mspminer-config",
	                      desc = "the configuration file used by mspminer",
	                      default = None)
	workflow.add_argument("bypass-clustering",
	                      desc = "do not cluster proteins into protein families",
	                      action = "store_true")
	workflow.add_argument("bypass-global-homology",
	                      desc = "do not annotate protein families based on global homology information",
	                      action = "store_true")
	workflow.add_argument("bypass-domain-motif",
	                      desc = "do not annotate protein families based on domain/motif information",
	                      action = "store_true")
	workflow.add_argument("bypass-interproscan",
	                      desc = "do not annotate protein families based on interproscan",
	                      action = "store_true")
	workflow.add_argument("bypass-pfam_to_go",
	                      desc = "do not annotate protein families based on pfam2go",
	                      action = "store_true")
	workflow.add_argument("bypass-domine",
	                      desc = "do not annotate protein families based on DOMINE database",
	                      action = "store_true")
	workflow.add_argument("bypass-sifts",
	                      desc = "do not annotate protein families based on SIFTS database",
	                      action = "store_true")
	workflow.add_argument("bypass-expatlas",
	                      desc = "do not annotate protein families based on Expression Atlas database",
	                      action = "store_true")
	workflow.add_argument("bypass-psortb",
	                      desc = "do not annotate protein families based on psortb",
	                      action = "store_true")
	workflow.add_argument("bypass-abundance",
	                      desc = "do not annotate protein families based on abundance information",
	                      action = "store_true")
	workflow.add_argument("split-number",
	                      desc="indicates number of spliting files for annotation based on sequence information",
	                      default = None)
	workflow.add_argument("bypass-integration",
	                      desc = "do not integrate annotations for protein families",
	                      action = "store_true")
	workflow.add_argument("output",
	                      desc = "provide an output folder which the workflow database and log is written. By default, thet be written to the users home directory",
	                      default = None)

	return workflow



def get_method_config (config_file):
	'''
	:param config_file:
	:return: config info for each analysis
	'''
	print('get_method_config')

	config_items = config.read_user_edit_config_file(config_file)
	family_conf = {}
	domain_motif_conf = {}
	abundance_conf = {}
	integration_conf = {}
	values = ["yes", "no", "Yes", "No"]

	if "global_homology" in config_items:
		for name in config_items["global_homology"].keys():
			myvalue = config_items["global_homology"][name]
			if not myvalue in values:
				print('The config value can not be recognized. Please check your config file!')
				continue
			family_conf[name] = myvalue
		# for each method

	if "domain_motif" in config_items:
		for name in config_items["domain_motif"].keys():
			myvalue = config_items["domain_motif"][name]
			if not myvalue in values:
				print('The config value can not be recognized. Please check your config file!')
				continue
			domain_motif_conf[name] = myvalue
		# for each method

	if "abundance" in config_items:
		for name in config_items["abundance"].keys():
			myvalue = config_items["abundance"][name]
			abundance_conf[name] = myvalue
		# for each method

	if "integration" in config_items:
		for name in config_items["integration"].keys():
			myvalue = config_items["integration"][name]
			if not myvalue in values:
				print('The config value can not be recognized. Please check your config file!')
				continue
			integration_conf[name] = myvalue
		# for each method

	return family_conf, domain_motif_conf, abundance_conf, integration_conf



def main(workflow):
	'''
	Characterization of protein families
	'''

	# get arguments
	args = workflow.parse_args()
	tmps_dir = os.path.join (os.getcwd(), "temp")
	tmps_dir = os.path.abspath(tmps_dir)
	if os.path.isdir(tmps_dir):
		os.system("rm -rf " + tmps_dir)

	# get configuration info
	default_characterization_conf = os.path.join(config.config_directory, "characterization.cfg")
	if args.characterization_config:
		args.characterization_config = os.path.abspath(args.characterization_config)
	else:
		args.characterization_config = default_characterization_conf
	print(args.characterization_config)
	
	# update configurations of characterization
	family_conf, domain_motif_conf, abundance_conf, integration_conf = get_method_config(args.characterization_config)
	if args.threads:
		args.threads = int(args.threads)
	else:
		args.threads = int(config.threads)
	if args.split_number:
		args.split_number = int(args.split_number)
	else:
		args.split_number = int(config.split_number)
	if args.mspminer_config:
		abundance_conf["mspminer"] = os.path.abspath(args.mspminer_config)
	else:
		abundance_conf["mspminer"] = config.mspminer
	if args.bypass_interproscan:
		domain_motif_conf["interproscan"] = "no"
		domain_motif_conf["pfam2go"] = "no"
		domain_motif_conf["domine"] = "no"
		domain_motif_conf["sifts"] = "no"
		domain_motif_conf["expatlas"] = "no"
	if args.bypass_pfamtogo:
		domain_motif_conf["pfam2go"] = "no"
	if args.bypass_domine:
		domain_motif_conf["domine"] = "no"
	if args.bypass_sifts:
		domain_motif_conf["sifts"] = "no"
	if args.bypass_expatlas:
		domain_motif_conf["expatlas"] = "no"
	if args.bypass_psortb:
		domain_motif_conf["psortb"] = "no"
	if "".join(config.phenotype) == "none":
		abundance_conf["dna_da"] = "no"


	# input and output folder
	#input_dir = args.input
	#input_dir = os.path.abspath(input_dir)
	output_dir = config.annotation_dir
	output_dir = os.path.abspath(output_dir)
	if not args.output:
		args.output = config.working_dir

	# get all input files
	#gene_catalog = config.gene_catalog
	gene_catalog_seq = config.gene_catalog_prot
	gene_catalog_count = config.gene_catalog_count

	# get all output files
	protein_family = config.protein_family
	protein_family_seq = config.protein_family_prot_seq
	protein_family_relab = config.protein_family_relab
	protein_family_ann = config.protein_family_ann
	protein_family_attr = config.protein_family_attr
	protein_family_ann_list = {}
	protein_ann_list = {}
	uniref_taxonomy_family = os.path.join(output_dir, config.basename + "_proteinfamilies_annotation.uniref90_annotation.tsv")
	uniref_taxonomy = os.path.join(output_dir, config.basename + "_protein_annotation.uniref90_annotation.tsv")
	taxonomy_annotation_family = os.path.join(output_dir, config.basename + "_proteinfamilies_annotation.taxonomy.tsv")
	taxonomy_annotation = os.path.join(output_dir, config.basename + "_protein_annotation.taxonomy.tsv")


	### STEP #1: clustering ###
	# if clustering action is provided, then cluster proteins into protein families
	if not args.bypass_clustering:
		print("Run clustering")
		myprotein_family, myprotein_family_output_folder = characterization.clustering (workflow, gene_catalog_seq,
		                                                                                args.threads,
		                                                                                output_dir, protein_family, protein_family_seq)

	### STEP #2: global-homology annotation ###
	# if global-homology action is provided, then directly extract annotations from database (UniProt)
	if not args.bypass_global_homology:
		print("Run global-homology annotation")
		myprotein_family_ann, myprotein_ann, homology_output_folder = characterization.global_homology_annotation (workflow, family_conf,
		                                                                                                          gene_catalog_seq,
		                                                                                                          args.threads,
		                                                                                                          output_dir, uniref_taxonomy_family, uniref_taxonomy,
		                                                                                                          protein_family_ann_list, protein_ann_list, protein_family, protein_family_seq)

	### STEP #3: domain-motif annotation ###
	# if domain-motif action is provided, then do annotations based on sequence information
	if not args.bypass_domain_motif:
		print("Run domain_motif annotation")
		myprotein_family_ann, myprotein_ann, sequence_output_folder = characterization.domain_motif_annotation (workflow, domain_motif_conf,
		                                                                                                          gene_catalog_seq,
		                                                                                                          args.split_number, args.threads,
		                                                                                                          output_dir, protein_family_ann_list, protein_ann_list, protein_family, protein_family_seq)


	### STEP #4: abundance annotation ###
	# if abundance action is provided, then do annotations based on abundance information
	if not args.bypass_abundance:
		print("Run abundance annotation")
		myprotein_family_ann, myprotein_ann, abundance_output_folder = characterization.abundance_annotation (workflow, abundance_conf,
		                                                                                                             gene_catalog_seq, gene_catalog_count,
		                                                                                                            uniref_taxonomy_family, uniref_taxonomy,
		                                                                                                            args.split_number, args.threads,
		                                                                                                            output_dir, protein_family, protein_family_relab, taxonomy_annotation_family, taxonomy_annotation,
		                                                                                                            protein_family_ann_list, protein_ann_list)

	### STEP #4: integrate annotations ###
	if not args.bypass_integration:
		print("Run integration annotation")
		myprotein_family_ann, myprotein_family_attr, annotation_output_folder = characterization.integration_annotation (workflow, integration_conf,
		                                                                                                               protein_family_ann_list, protein_ann_list,
		                                                                                                               uniref_taxonomy_family, uniref_taxonomy,
		                                                                                                               taxonomy_annotation_family, taxonomy_annotation,
		                                                                                                               args.threads,
		                                                                                                               output_dir, protein_family_ann, protein_family_attr)

	### start the workflow
	workflow.go()



if __name__ == "__main__":
	main(parse_cli_arguments())




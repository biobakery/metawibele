#!/usr/bin/env python

"""
Preprocessing workflow:
1) assembly metagenomic shotgun sequencing reads
2) identify ORFs
3) build non-redundancy gene catalogs
4) calculate reads counts table for gene catalogs

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

import fnmatch
import logging
import os
import sys

# import the workflow class from anadama2
from anadama2 import Workflow

# import the utilities functions and config settings from MetaWIBELE
try:
	from metawibele import utilities, config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
		         " Please check your install.")

# import the library of MetaWIBELE tasks for characterization
try:
	from metawibele.tools import preprocessing_tasks
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
		         " Please check your install.")

VERSION = config.version

def parse_cli_arguments():
	'''
	 Parses any command-line arguments passed into the workflow.
	the version number will appear when running this script with the "--version" option
	the description will appear when running this script with the "--help" option
	create a workflow instance, providing the version number and description
	'''

	workflow = Workflow(version=VERSION, description="A workflow to preprocess shotgun sequencing reads of metagenomes " 
													"with tasks of metagenomic assembly, gene calling, " 
												"building gene catalogs and generating gene abundance for each sample.")

	# add the custom arguments to the workflow
	workflow.add_argument("threads",
	                      desc = "number of threads/cores for each task to use",
	                      default = None)
	workflow.add_argument("extension-paired",
	                      desc = "provide the extension for paired fastq files using comma to separate, e.g. .R1.fastq.gz,.R2.fastq.gz | .R1.fastq,.R2.fastq",
						  default = None) 
	workflow.add_argument("extension",
	                      desc = "provide the extension for all fastq files",
						  choices = [".fastq.gz", ".fastq"],
	                      default = ".fastq.gz")
	workflow.add_argument("gene-call-type",
	                      desc = "specify which type of gene calls will be used",
						  choices = ['prokka', 'prodigal', 'both'],
	                      default = 'both')
	workflow.add_argument("bypass-assembly",
	                      desc = "do not run assembly",
	                      action = "store_true")
	workflow.add_argument("bypass-gene-calling",
	                      desc = "do not call ORFs",
	                      action = "store_true")
	workflow.add_argument("bypass-gene-catalog",
	                      desc = "do not build gene catalogs",
	                      action = "store_true")
	workflow.add_argument("output-basename",
	                      desc = "provide the basename for output files",
	                      default = "mgx")

	return workflow


def main(workflow):
	'''
	Characterization of protein families
	'''

	# get arguments
	args = workflow.parse_args()
	if args.threads:
		args.threads = int(args.threads)
	else:
		args.threads = int(config.threads)
	if args.output_basename:
		args.output_basename = args.output_basename
	else:
		args.output_basename = config.basename

	# input and output folder
	input_dir = args.input # reads fastq files
	output_dir = args.output
	input_dir = os.path.abspath(input_dir)
	output_dir = os.path.abspath(output_dir)
	extension_paired = args.extension_paired
	file_extension = args.extension

	# get all output files
	assembly_dir = os.path.join(output_dir, "assembly/")
	contigs = os.path.join(output_dir, args.output_basename + "_contig_sequence.fasta")
	assembly_extentsion = ".contigs.fa"

	prokka_dir = os.path.join(output_dir, "gene_annotation/")
	prodigal_dir = os.path.join(output_dir, "gene_calls/")
	gene_file = os.path.join(output_dir, args.output_basename + "_combined_gene.fna")
	gene_PC_file = os.path.join(output_dir, args.output_basename + "_combined_gene_protein_coding.sorted.fna")
	protein_file = os.path.join(output_dir, args.output_basename + "_combined_protein.faa")
	protein_sort = os.path.join(output_dir, args.output_basename + "_combined_protein.sorted.faa")
	gene_info = os.path.join(output_dir, args.output_basename + "_gene_info.tsv")
	complete_gene = os.path.join(output_dir, args.output_basename + "_combined_gene_protein_coding.complete.sorted.fna")
	complete_protein = os.path.join(output_dir, args.output_basename + "_combined_protein.complete.sorted.faa")

	mapping_dir = os.path.join(output_dir, "mapping/")
	prefix_gene_catalog = os.path.join(output_dir, args.output_basename + "_genecatalogs.centroid")
	gene_catalog = os.path.join(output_dir, args.output_basename + "_genecatalogs.clstr")
	gene_catalog_nuc = os.path.join(output_dir, args.output_basename + "_genecatalogs.centroid.fna")
	gene_catalog_prot = os.path.join(output_dir, args.output_basename + "_genecatalogs.centroid.faa")
	gene_catalog_saf = os.path.join(output_dir, args.output_basename + "_genecatalogs.centroid.saf.gtf")
	gene_catalog_count = os.path.join(output_dir, args.output_basename + "_genecatalogs_counts.all.tsv")

	
	### STEP #1: assembly ###
	# if assembly action is provided, then assembly MGX reads
	if not args.bypass_assembly:
		mycontigs = preprocessing_tasks.assembly (workflow, input_dir,
		                                         args.extension, extension_paired,
		                                         args.threads,
		                                         assembly_dir, contigs)

	### STEP #2: gene calling ###
	# if gene-calling action is provided, then identify ORFs
	if not args.bypass_gene_calling:
		mygene, myprotein = preprocessing_tasks.gene_calling (workflow, assembly_dir, assembly_extentsion,
															 input_dir, args.extension, extension_paired,
		                                                     args.gene_call_type, prokka_dir, prodigal_dir,
		                                                     args.threads,
		                                                     gene_file, gene_PC_file, protein_file, protein_sort,
		                                                     gene_info, complete_gene, complete_protein)

	### STEP #3: gene-catalog  ###
	# if gene-catalog action is provided, then build gene catalogs and calculate abundance per gene catalogs across samples
	if not args.bypass_gene_catalog:
		mygene_catalog, mycounts = preprocessing_tasks.gene_catalog (workflow, 
																	 complete_gene, complete_protein,
																	 input_dir, args.extension, extension_paired,
		                                                             args.threads,
		                                                             prefix_gene_catalog, gene_catalog,
		                                                             gene_catalog_nuc, gene_catalog_prot,
		                                                             mapping_dir, gene_catalog_saf, gene_catalog_count)

	### start the workflow
	workflow.go()

if __name__ == "__main__":
	main(parse_cli_arguments())

#!/usr/bin/env python3

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
from metawibele import utilities, config
# import the library of MetaWIBELE tasks for characterization
from metawibele.tools import preprocessing_tasks

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
	                      desc="number of threads/cores for each task to use",
	                      default=20)
	workflow.add_argument("sample-file",
	                      desc=" the list file of samples",
	                      required=True)
	workflow.add_argument("extension-paired",
	                      desc="indicates the extension for paired fastq files, e.g. _1.fastq.gz,_2.fastq.gz",
	                      required=True)
	workflow.add_argument("extension-orphan",
	                      desc="indicates the extension for orphan fastq files",
	                      default="none")
	workflow.add_argument("assembly",
	                      desc="indicates whether or not do assembly",
	                      default=True)
	workflow.add_argument("gene-calling",
	                      desc="indicates whether or not call ORFs",
	                      default=True)
	workflow.add_argument("gene-catalog",
	                      desc="indicates whether or not build gene catalogs",
	                      default=True)
	workflow.add_argument("output-basename",
	                      desc="specify the basename for output files",
	                      default="metawibele")

	return workflow


def main(workflow):
	'''
	Characterization of protein families
	'''

	# get arguments
	args = workflow.parse_args()

	# input and output folder
	input_dir = args.input # reads fastq files
	output_dir = args.output
	extension_paired = args.extension_paired
	extension_orphan = args.extension_orphan
	file_extension = "none"
	if extension_paired != "none":
		file_extension = extension_paired
	if extension_orphan != "none":
		if file_extension != "none":
			file_extension = file_extension + "," + extension_orphan
		else:
			file_extension = extension_orphan

	# get all output files
	assembly_dir = output_dir + "/assembly/"
	contigs = output_dir + "/" + args.basename + "_contig_sequence.fasta"
	assembly_extentsion = ".contigs.fa"

	prokka_dir = output_dir + "/gene_annotation/"
	prodigal_dir = output_dir + "/gene_calls/"
	gene_file = output_dir + "/" + args.basename + "_combined_gene.fna"
	gene_PC_file = output_dir + "/" + args.basename + "_combined_gene_protein_coding.sorted.fna"
	protein_file = output_dir + "/" + args.basename + "_combined_protein.faa"
	protein_sort = output_dir + "/" + args.basename + "_combined_protein.sorted.faa"
	gene_info = output_dir + "/" + args.basename + "_gene_info.tsv"
	complete_gene = output_dir + "/" + args.basename + "_combined_gene_protein_coding.complete.sorted.fna"
	complete_protein = output_dir + "/" + args.basename + "_combined_protein.complete.sorted.faa"

	mapping_dir = output_dir + "/mapping/"
	prefix_gene_catalog = output_dir + "/" + args.basename + "_genecatalogs.centroid"
	gene_catalog = output_dir + "/" + args.basename + "_genecatalogs.clstr"
	gene_catalog_nuc = output_dir + "/" + args.basename + "_genecatalogs.centroid.fna"
	gene_catalog_prot = output_dir + "/" + args.basename + "_genecatalogs.centroid.faa"
	gene_catalog_saf = output_dir + "/" + args.basename + "_genecatalogs.centroid.saf.gtf"
	gene_catalog_count = output_dir + "/" + args.basename + "_genecatalogs_counts.all.tsv"

	
	### STEP #1: assembly ###
	# if assembly action is provided, then assembly MGX reads
	if args.assembly == True:
		mycontigs = preprocessing_tasks.assembly (workflow, input_dir, args.sample_file,
		                                         extension_paired, extension_orphan,
		                                         args.threads,
		                                         assembly_dir, contigs)

	### STEP #2: gene calling ###
	# if gene-calling action is provided, then identify ORFs
	if args.gene_calling == True:
		mygene, myprotein = preprocessing_tasks.gene_calling (workflow, assembly_dir, assembly_extentsion,
		                                                     args.sample_file,
		                                                     prokka_dir, prodigal_dir,
		                                                     args.threads,
		                                                     gene_file, gene_PC_file, protein_file, protein_sort,
		                                                     gene_info, complete_gene, complete_protein)

	### STEP #3: gene-catalog  ###
	# if gene-catalog action is provided, then build gene catalogs and calculate abundance per gene catalogs across samples
	if args.gene_catalog == True:
		mygene_catalog, mycounts = preprocessing_tasks.gene_catalog (workflow, complete_gene, complete_protein,
		                                                             input_dir, args.sample_file, file_extension,
		                                                             args.threads,
		                                                             prefix_gene_catalog, gene_catalog,
		                                                             gene_catalog_nuc, gene_catalog_prot,
		                                                             mapping_dir, gene_catalog_saf, gene_catalog_count)

	### start the workflow
	workflow.go()

if __name__ == "__main__":
	main(parse_cli_arguments())

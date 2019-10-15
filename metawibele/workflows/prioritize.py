#!/usr/bin/env python

"""
MeteWIBELE workflow: MeteWIBELE prioritization workflow
1) unsupervised prioritization based on abundance and prevalence
2) supervised prioritization based on association with phenotypes
3) selected prioritization based on sequence annotations
4) format and finalize prioritization


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

# import the library of MetaWIBELE tasks for prioritization
from metawibele.tasks import prioritization

# import the utilities functions and config settings from MetaWIBELE
from metawibele import utilities, config


def parse_cli_arguments ():
	'''
	 Parses any command-line arguments passed into the workflow.
	the version number will appear when running this script with the "--version" option
	the description will appear when running this script with the "--help" option
	create a workflow instance, providing the version number and description
	'''

	workflow = Workflow(version = "0.0.1", description = "A workflow for MetaWIBELE prioritization")

	# add the custom arguments to the workflow
	workflow.add_argument("threads",
	                      desc = "number of threads/cores for each task to use",
	                      default = 20)
	workflow.add_argument("quantification-based",
	                      desc = "indicates whether or not prioritize protein families based on quantitative criteria",
	                      default = True,
	                      action = "store_true")
	workflow.add_argument("quantitative-config",
	                      desc = "the configuration file of quantitative prioritization",
	                      default = "prioritization.cfg"),
	workflow.add_argument("annotation-based",
	                      desc = "indicates whether or not prioritize protein families based on sequence annotation",
	                      default = True,
	                      action = "store_true")
	workflow.add_argument("selection-config",
	                      desc = "the configuration file of annotation-based prioritization",
	                      default = "selection.cfg"),
	workflow.add_argument("selection-file",
	                      desc="the output file name of annotation-based prioritization",
	                      default = "selected_prioritization.tsv")
	workflow.add_argument("finalized",
	                      desc = "indicates whether or not finalize prioritized protein families",
	                      default = True,
	                      action = "store_true")

	return workflow


def main(workflow):
	'''
	Prioritization of protein families
	'''

	# get arguments
	args = workflow.parse_args()

	# input and output folder
	input_dir = args.input
	output_dir = config.priority_dir

	# collect input files
	protein_family = config.protein_family
	protein_family_seq = config.protein_family_seq
	protein_family_relab = config.protein_family_relab
	protein_family_ann = config.protein_family_ann
	protein_family_attr = config.protein_family_attr

	# output files
	unsupervised_rank = priority_dir + "/" + basename + "_unsupervised_prioritization.rank.tsv"
	supervised_rank = priority_dir + "/" + basename + "_supervised_prioritization.rank.tsv"
	unsupervised_priority = priority_dir + "/" + basename + "_unsupervised_prioritization.priority.tsv"
	supervised_priority = priority_dir + "/" + basename + "_supervised_prioritization.priority.tsv"
	selected_priority = priority_dir + "/" + args.selection_file

	final_unsupervised_rank = priority_dir + "/" + basename + "_unsupervised_prioritization.rank.table.tsv"
	final_supervised_rank = priority_dir + "/" + basename + "_supervised_prioritization.rank.table.tsv"
	final_unsupervised_priority = priority_dir + "/" + basename + "_unsupervised_prioritization.priority.table.tsv"
	final_supervised_priority = priority_dir + "/" + basename + "_supervised_prioritization.priority.table.tsv"
	final_selected_priority = re.sub(".tsv", ".table.tsv", selected_priority)


	### STEP #1: quantification-based prioritization ###
	# if quantification-based action is provided, then prioritize protein families using quantitative criteria
	if args.quantification_based:
		unsupervised_file, supervised_file = prioritization.quantification_based_prioritization (workflow, args.quantitative_config,
		                                                                        protein_family_ann, protein_family_attr,
		                                                                        output_dir)


	### STEP #2: annotation-based prioritization ###
	# if annotation-based action is provided, then prioritize protein families based on sequence annotations (selection factor)
	if args.annotation_based:
		myselection = prioritization.annotation_based_prioritization (workflow, args.selection_config,
		                                                             protein_family_ann,
		                                                             supervised_priority,
		                                                             output_dir, selected_priority)

	### STEP #3: finalized annotation ###
	# if finalized action is provided, then format and fianlize prioritizations
	if args.finalized:
		prioritization.finalize_prioritization (workflow,
		                                        unsupervised_rank, unsupervised_priority,
		                                        supervised_rank, supervised_priority, selected_priority,
		                                        output_dir,
		                                        final_unsupervised_rank, final_unsupervised_priority,
		                                        final_supervised_rank, final_supervised_priority, final_selected_priority)

	### start the workflow
	workflow.go()



if __name__ == "__main__":
	main(parse_cli_arguments())


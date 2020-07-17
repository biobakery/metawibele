#!/usr/bin/env python

"""
MeteWIBELE workflow: MeteWIBELE prioritization workflow
1) unsupervised prioritization based on abundance and prevalence
2) supervised prioritization based on association with phenotypes
3) selected prioritization based on interested functions
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
import re

# import the workflow class from anadama2
from anadama2 import Workflow

# import the library of MetaWIBELE tasks for prioritization
try:
	from metawibele.tasks import prioritization
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

	workflow = Workflow(version = VERSION, description = "A workflow for MetaWIBELE prioritization", remove_options=["input","output"])

	# add the custom arguments to the workflow
	workflow.add_argument("threads",
						desc = "number of threads/cores for each task to use",
						default = None)
	workflow.add_argument("prioritization-config",
	                    desc = "the configuration file for prioritization",
	                    default = None)
	workflow.add_argument("vignette-config",
	                    desc = "the file with specific functions of interest used as binary filtering for prioritization",
	                    default = "none")
	workflow.add_argument("bypass-mandatory",
	                     desc = "do not prioritize protein families based on quantitative criteria (mandatory prioritization)",
						 action = "store_true")
	workflow.add_argument("bypass-optional",
	                     desc = "do not prioritize protein families based on selecting our for interested annotations (optional prioritization)",
	                     action = "store_true")
	workflow.add_argument("bypass-finalized",
	                     desc = "do not finalize prioritized protein families",
						 action = "store_true")
	workflow.add_argument("selected-output",
	                    desc = "the output file name for the prioritized protein families by binary filtering",
	                    default = None)
	workflow.add_argument("output",
	                    desc = "provide an output folder which the workflow database and log is written. By default, thet be written to the users home directory",
	                    default = None)

	return workflow


def main(workflow):
	'''
	Prioritization of protein families
	'''

	# get arguments
	args = workflow.parse_args()
	if args.threads:
		args.threads = int(args.threads)
	else:
		args.threads = int(config.threads)

	# input and output folder
	#input_dir = args.input
	#input_dir = os.path.abspath(input_dir)
	priority_dir = config.priority_dir
	priority_dir = os.path.abspath(priority_dir)
	if not args.output:
		args.output = config.working_dir 


	# get config file
	default_prioritization_conf = os.path.join(config.config_directory, "prioritization.cfg")
	if args.prioritization_config:
		args.prioritization_config = os.path.abspath(args.prioritization_config)
	else:
		args.prioritization_config = default_prioritization_conf
	print(args.prioritization_config)

	# collect input files
	protein_family = config.protein_family
	protein_family_seq = config.protein_family_prot_seq
	protein_family_relab = config.protein_family_relab
	protein_family_ann = config.protein_family_ann
	protein_family_attr = config.protein_family_attr

	# output files
	unsupervised_rank = os.path.join(priority_dir, config.basename + "_unsupervised_prioritization.rank.tsv")
	supervised_rank = os.path.join(priority_dir, config.basename + "_supervised_prioritization.rank.tsv")
	if args.selected_output:
		selected_priority = os.path.join(priority_dir, os.path.basename(args.selected_output))
	else:	
		selected_priority = os.path.join(priority_dir, config.basename + "_supervised_prioritization.rank.selected.tsv") 

	final_unsupervised_rank = os.path.join(priority_dir, config.basename + "_unsupervised_prioritization.rank.table.tsv")
	final_supervised_rank = os.path.join(priority_dir, config.basename + "_supervised_prioritization.rank.table.tsv")
	final_selected_priority = re.sub(".tsv", ".table.tsv", selected_priority)


	### STEP #1: mandatory prioritization: quantification-based ranking ###
	# if mandatory action is provided, then prioritize protein families using quantitative criteria
	if not args.bypass_mandatory:
		unsupervised_file, supervised_file = prioritization.mandatory_prioritization (workflow, args.prioritization_config,
		                                                                        protein_family_ann, protein_family_attr,
		                                                                        priority_dir)


	### STEP #2: optional prioritization: binary filtering ###
	# if optional action is provided, then prioritize protein families based on interested functions (selection factor)
	if not args.bypass_optional:
		if not "".join(config.phenotype) == "none":
			myselection = prioritization.optional_prioritization (workflow, args.prioritization_config, args.vignette_config,
		                                                             protein_family_ann,
		                                                             supervised_rank,
		                                                             priority_dir, selected_priority)
	

	### STEP #4: finalized annotation ###
	# if finalized action is provided, then format and fianlize prioritizations
	if not args.bypass_finalized:
		prioritization.finalize_prioritization (workflow,
		                                        unsupervised_rank,
		                                        supervised_rank, selected_priority,
		                                        priority_dir,
		                                        final_unsupervised_rank,
		                                        final_supervised_rank, final_selected_priority)

	### start the workflow
	workflow.go()



if __name__ == "__main__":
	main(parse_cli_arguments())



#!/usr/bin/env python

"""
MetaWIBELE Workflow: prioritization module
A collection of tasks for workflows with prioritization

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
import subprocess
import itertools

from anadama2.tracked import TrackedExecutable, TrackedDirectory

# import the utilities functions and config settings from MetaWIBELE
from metawibele import utilities, config, files


def mandatory_prioritization (workflow, prioritization_conf,
                                         protein_family_ann, protein_family_attr,
                                         output_folder):
	"""
	This set of tasks will run prioritization using quantitative criteria.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		prioritization_conf: Configuration file for quantitative prioritization.
		protein_family_ann: Finalized annotation file for protein .
		protein_family_attr: Finalized attribue file for annotations.

	Requires:
		config file
		annotation files

	Returns:
		string: the name of prioritized file.

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add quantification_based_prioritization tasks
		myrank, mypriority = prioritization.mandatory_prioritization (workflow, args.prioritization_conf,
		                                                                        protein_family_ann, protein_family_attr,
		                                                                        output_dir)
		# run the workflow
		workflow.go()
	"""

	# get the clustering output files
	priority_dir = output_folder
	unsupervised_rank = priority_dir + "/" + basename + "_unsupervised_prioritization.rank.tsv"
	supervised_rank = priority_dir + "/" + basename + "_supervised_prioritization.rank.tsv"
	unsupervised_priority = priority_dir + "/" + basename + "_unsupervised_prioritization.priority.tsv"
	supervised_priority = priority_dir + "/" + basename + "_supervised_prioritization.priority.tsv"


	# run unsupervised prioritization
	mylog = re.sub(".tsv", ".log", unsupervised_rank)
	workflow.add_task(
			"quantify_prioritization.py -c [depends[0]] -m unsupervised -w fixed -a [depends[1]] -b [depends[2]] -o [args[0]] > [args[1]] 2>&1",
			depends = [prioritization_conf, protein_family_ann, protein_family_attr, TrackedExecutable("quantify_prioritization.py")],
			targets = [unsupervised_rank, unsupervised_priority],
			args = [priority_dir, mylog],
			name = "quantify_prioritization")

	# run supervised prioritization
	mylog = re.sub(".tsv", ".log", supervised_rank)
	workflow.add_task(
			"quantify_prioritization.py -c [depends[0]] -m supervised -w fixed -a [depends[1]] -b [depends[2]] -o [args[0]] > [args[1]] 2>&1",
			depends = [prioritization_conf, protein_family_ann, protein_family_attr, TrackedExecutable("quantify_prioritization.py")],
			targets = [supervised_rank, supervised_priority],
			args = [priority_dir, mylog],
			name = "quantify_prioritization")

	return unsupervised_priority, supervised_priority


def optional_prioritization (workflow, prioritization_conf,
                                     protein_family_ann, supervised_priority,
                                     output_folder, selected_priority):
	"""
	This set of tasks will run prioritization using functional annotations as optional fiters.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		prioritization_conf: Configuration file for quantitative prioritization.
		protein_family_ann: Finalized annotation file for protein .
		supervised_priority: Supervised prioritization file.

	Requires:
		config file
		annotation and quantitative prioritization files

	Returns:
		string: the name of annotation-based prioritization file.

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# annotation_based_prioritization tasks
		myselection = prioritization.optional_prioritization (workflow, args.prioritization_conf,
		                                                             protein_family_ann,
		                                                             supervised_priority,
		                                                             output_dir, selected_priority)
		# run the workflow
		workflow.go()
	"""

	# get the clustering output files
	priority_dir = output_folder

	# run annotation-based prioritization
	mylog = re.sub(".tsv", ".log", selected_priority)
	workflow.add_task(
			"filter_prioritization.py -c [depends[0]] -a [depends[1]] -p [depends[2]] -o [targets[0]] > [args[0]] 2>&1",
			depends = [prioritization_conf, protein_family_ann, supervised_priority, TrackedExecutable("filter_prioritization.py")],
			targets = [selected_priority],
			args = [mylog],
			name = "filter_prioritization")

	return selected_priority


def finalize_prioritization (workflow,
                             unsupervised_rank, unsupervised_priority,
                             supervised_rank, supervised_priority, selected_priority,
                             output_folder,
                             final_unsupervised_rank, final_unsupervised_priority,
                             final_supervised_rank, final_supervised_priority, final_selected_priority):

	"""
	This set of tasks will format prioritization files

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		raw prioritized files
		finalized prioritized files

	Requires:
		raw prioritized files


	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add quality control tasks for the fastq files
		finalize_prioritization (workflow,
                             unsupervised_rank, unsupervised_priority,
                             supervised_rank, supervised_priority, selected_priority,
                             output_folder,
                             final_unsupervised_rank, final_unsupervised_priority,
                             final_supervised_rank, final_supervised_priority, final_selected_priority)
		# run the workflow
		workflow.go()
	"""

	# format prioritization
	mylog = re.sub(".tsv", ".log", final_unsupervised_rank)
	workflow.add_task(
			"finalize_prioritization.py -i [depends[0]] -o [targets[0]] > [args[0]] 2>&1",
			depends = [unsupervised_rank, TrackedExecutable("finalize_prioritization.py")],
			targets = [final_unsupervised_rank],
			args = [mylog],
			name = "finalize_prioritization")

	mylog = re.sub(".tsv", ".log", final_unsupervised_priority)
	workflow.add_task(
			"finalize_prioritization.py -i [depends[0]] -o [targets[0]] > [args[0]] 2>&1",
			depends = [unsupervised_priority, TrackedExecutable("finalize_prioritization.py")],
			targets = [final_unsupervised_priority],
			args = [mylog],
			name = "finalize_prioritization")

	mylog = re.sub(".tsv", ".log", final_supervised_rank)
	workflow.add_task(
			"finalize_prioritization.py -i [depends[0]] -o [targets[0]] > [args[0]] 2>&1",
			depends = [supervised_rank, TrackedExecutable("finalize_prioritization.py")],
			targets = [final_supervised_rank],
			args = [mylog],
			name = "finalize_prioritization")

	mylog = re.sub(".tsv", ".log", final_supervised_priority)
	workflow.add_task(
			"finalize_prioritization.py -i [depends[0]] -o [targets[0]] > [args[0]] 2>&1",
			depends = [supervised_priority, TrackedExecutable("finalize_prioritization.py")],
			targets = [final_supervised_priority],
			args = [mylog],
			name = "finalize_prioritization")

	mylog = re.sub(".tsv", ".log", final_selected_priority)
	workflow.add_task(
			"finalize_prioritization.py -i [depends[0]] -o [targets[0]] > [args[0]] 2>&1",
			depends = [selected_priority, TrackedExecutable("finalize_prioritization.py")],
			targets = [final_selected_priority],
			args = [mylog],
			name = "finalize_prioritization")


def moduled_prioritization (workflow, prioritization_conf,
                            protein_family_ann, supervised_priority,
							output_folder, moduled_priority):

	"""
	This set of tasks will identify potential functional modules for prioritized proteins

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		prioritization_conf: Configuration file for quantitative prioritization.
		protein_family_ann: Finalized annotation file for protein .
		supervised_priority: Supervised prioritization file.

	Requires:
		config file
		annotation and quantitative prioritization files

	Returns:
		string: the name of annotation-based prioritization file.

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# annotation_based_prioritization tasks
		mymodule = prioritization.moduled_prioritization (workflow, args.prioritization_conf,
		                                                             protein_family_ann,
		                                                             supervised_priority,
		                                                             output_dir, moduled_priority)

		workflow.go()
	"""

	# get the clustering output files
	priority_dir = output_folder

	# run moduling prioritization
	mylog = re.sub(".tsv", ".log", moduled_priority)
	workflow.add_task(
			"cluster_prioritization.py -c [depends[0]] -a [depends[1]] -p [depends[2]] -o [targets[0]] > [args[0]] 2>&1",
			depends = [prioritization_conf, protein_family_ann, supervised_priority, TrackedExecutable("cluster_prioritization.py")],
			targets = [selected_priority],
			args = [mylog],
			name = "cluster_prioritization")

	return moduled_priority


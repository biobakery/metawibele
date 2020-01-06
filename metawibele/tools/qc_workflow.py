#!/usr/bin/env python

"""
mgx_QC_workflow.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run quality control of raw data

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

import itertools
import os
import re
import shutil
import tempfile
from anadama2 import Workflow
from glob2 import glob
from itertools import chain
from anadama2.tracked import TrackedExecutable, TrackedDirectory

# import the utilities functions and config settings from MetaWIBELE
from metawibele import utilities, config

VERSION = config.version

def parse_cli_arguments():
	"""Run kneaddata
		This set of tasks will run kneaddata on the input files provided. It will run with single-end or paired-end input files.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		input (string): The path to fastq files for input to kneaddata.
        sample_file (string): A file for sample names.
		output (string): The path of the output folder.
		trimmomatic_options (string): trimmomatic options when running kneaddata (optional).
		additional_options (string): Additional options when running kneaddata (optional).
		remove_intermediate_output (bool): Remove intermediate output files.
		contaminant_db (string): The databases to use with kneaddata (optional).
								Allow for a single path or multiple paths in one string comma-delimited.
		file_extension (string): The extension for fastq files. 
		threads (int): The number of threads/cores for kneaddata to use.

	Requires:
		kneaddata v0.7.0+: A tool to perform quality control on metagenomic and metatranscriptomic sequencing data
	
	Returns:
		None
	
	Example:
    	python qc_workflow.py --input /my/reads/path --sample-file sample.txt --output myoutput --file-extension "_R1.fastq.gz,_R2.fastq.gz"
	
	"""
	
	workflow = Workflow(version=VERSION, description='A workflow to run kneaddata on the input files provided '
	                                               'to perform quality control.')
	workflow.add_argument('sample-file', 
						desc='Sample files including sample names (string).',
						required=True)
	workflow.add_argument("trimmomatic-options",
						desc="options for trimmomatic (string): ILLUMINACLIP:/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 (optional)",
	                    default="none")
	workflow.add_argument("additional-options",
						desc="additional_options (string): Additional options when running kneaddata (optional)",
						default="none")
	workflow.add_argument("remove-intermediate-output",
						desc="remove_intermediate_output (bool): Remove intermediate output files.",
						default=True)
	workflow.add_argument('contaminant-db',
	                      desc='Select reference sequences for the contamination you are trying to remove. '
	                           'It is KneadData databases including the indexed redernece sequences.',
	                      default="none")
	workflow.add_argument('file-extension', 
						desc='Extension of input fastq files (string)', 
						default=".R1.fastq.gz,.R2.fastq.gz")
	workflow.add_argument('threads', 
						desc='number of threads/cores for each task to use', 
						default=6)

	return workflow


def main(workflow):
	args = workflow.parse_args()
	threads = args.threads

	# collect samples names
	samples = []
	open_file = open(args.sample_file, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		samples.append(info[0])
	# foreach sample
	open_file.close()
	split_dir = args.input

	# get the kneaddata final output dir
	output_dir = os.path.join(args.output, "kneaddata")

	# collect input files
	paired = "flase"
	extension = args.file_extension.split(",")
	if len(extension) == 2:
		paired = "True"
	else:
		paired = "Flase"
	split_files = []
	for sample in samples:
		if paired == "True":
			f_seq = os.path.join(split_dir, sample + extension[0])
			r_seq = os.path.join(split_dir, sample + extension[1])
			f_seq_cleaned = os.path.join(output_dir, '%s.repeats.removed.1.fastq' % sample)
			r_seq_cleaned = os.path.join(output_dir, '%s.repeats.removed.2.fastq' % sample)
			f_seq_unmatched = os.path.join(output_dir, '%s.repeats.removed.unmatched.1.fastq' % sample)
			r_seq_unmatched = os.path.join(output_dir, '%s.repeats.removed.unmatched.2.fastq' % sample)
			# reorder the input files so they are a set of paired files
			input_files = [f_seq, r_seq]
			# output files
			output_log = os.path.join(output_dir, "%s_kneaddata.log" %sample)
			output_files = [f_seq_cleaned, r_seq_cleaned, f_seq_unmatched, r_seq_unmatched]
			kneaddata_output_repeats_removed_fastq = utilities.name_files(sample, output_dir, extension="repeats.removed.fastq")
			# add the second input file to the kneaddata arguments
			second_input_option=" --input [depends[1]] "
			# determine time/memory equations based on the two input files
			time_equation="3*6*60 if ( file_size('[depends[0]]') + file_size('[depends[1]]') ) < 10 else 5*6*60"
			mem_equation="3*12*1024 if ( file_size('[depends[0]]') + file_size('[depends[1]]') ) < 10 else 6*12*1024"
			rename_final_output = ""	
		else:
			f_seq = os.path.join(split_dir, sample + extension[0])
			f_seq_cleaned = os.path.join(output_dir, '%s.repeats.removed.fastq' % sample)
			# input file
			input_files = [f_seq]
			# output file
			output_log = os.path.join(output_dir, "%s_kneaddata.log" %sample)
			output_files = [f_seq_cleaned]
			kneaddata_output_repeats_removed_fastq = utilities.name_files(sample, output_dir, extension="repeats.removed.fastq")
			# the second input option is not used since these are single-end input files
			second_input_option=" "
			# determine time/memory equations based on the single input file
			time_equation="3*6*60 if file_size('[depends[0]]') < 10 else 5*6*60"
			mem_equation="3*12*1024 if file_size('[depends[0]]') < 10 else 6*12*1024"
			# need to rename the final output file here to the sample name
			rename_final_output = " && mv [args[3]] [targets[0]]"
		
		outstr = " ".join(output_files)
		split_files.append((sample, input_files, output_log, outstr, output_files, kneaddata_output_repeats_removed_fastq))

	
	# set additional options to empty string if not provided
	if args.additional_options == "none":
		additional_options = ""
	else:
		additional_options = args.additional_options
	# always run with the serial option, which presents read counts in log in the manner expected
	# by the visualization workflows (in serial filtering in the order of the databases provided)
	# also run tandem repeat filtering
	additional_options+=" --serial --run-trf "
    
	# add option to remove intermediate output, if set
	if args.remove_intermediate_output == True:
		additional_options += " --remove-intermediate-output "

	if args.trimmomatic_options == "none":
		trimmomatic_options = ""
	else:
		trimmomatic_options = " --trimmomatic-options  "  + args.trimmomatic_options
	
    # create the database command option string to provide zero or more databases to kneaddata
	if args.contaminant_db == "none":
		optional_arguments = ""
	elif isinstance(args.contaminant_db, list):
        # start the string with the kneaddata option and add an option for each database
		optional_arguments = " --reference-db "+" --reference-db ".join(args.contaminant_db)
	elif isinstance(args.contaminant_db, str) and "," in args.contaminant_db:
        # split the paths by comma
		database_list=list(filter(lambda x: x, args.contaminant_db.split(",")))
        # start the string with the kneaddata option and add an option for each database
		optional_arguments = " --reference-db "+" --reference-db ".join(database_list)        
	else:
		optional_arguments=" --reference-db " + args.contaminant_db
        
	
	# gzip fastq files
	rename_final_output = rename_final_output + " && gzip [args[4]]"

	# create a task for each set of input and output files to run kneaddata
	# rename file with repeats in name to only sample name
	os.system("mkdir -p " + output_dir)
	for (sample, depends, output_log, target_str, targets, intermediate_file) in split_files:
		gzipout = []
		for item in targets:
			myitem = re.sub(".fastq", ".fastq.gz", item)
			gzipout.append(myitem)
		seq_base = sample
		workflow.add_task_gridable(
            "kneaddata --input [depends[0]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]] " + second_input_option + optional_arguments + " " + additional_options + " " + trimmomatic_options + " > [args[5]] 2>&1 " + rename_final_output + " >> [args[5]] 2>&1 ",
			depends = utilities.add_to_list(depends, TrackedExecutable("kneaddata")),
			targets = gzipout,
			args = [output_dir, threads, sample, intermediate_file, target_str, output_log],
			time = time_equation, # 6 hours or more depending on file size
			mem = mem_equation, # 12 GB or more depending on file size
			cores = threads, # time/mem based on 8 cores
			name = utilities.name_task(sample, "kneaddata")) # name task based on sample name


	workflow.go()


if __name__ == "__main__":
	main(parse_cli_arguments())

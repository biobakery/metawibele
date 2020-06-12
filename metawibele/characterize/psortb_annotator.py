#!/usr/bin/env python

"""
psortb_annotator.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run PSORTb to predict secreted peptides based on representatives of clusters

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
import shutil
import tempfile
import re
import argparse

description = """
A workflow to annotate subcellular annotation
"""

def parse_cli_arguments():
	"""Parses any command-line arguments passed into the workflow.
	"""

	parser = argparse.ArgumentParser( description=description )
	parser.add_argument("--split-file", "-s",
						help = 'split files name',
						required = True)
	parser.add_argument("--threads", "-t",
						help = 'number of threads/cores for each task to use',
						default = 4,
						required = True)
	parser.add_argument("--input", "-i",
						help = 'input file',
						required = True)
	parser.add_argument("--output", "-o",
						help = 'output direcory',
						required = True)
	args = parser.parse_args()

	return args


def main():
	args = parse_cli_arguments()

	# ================================================
	# collect sequences
	# ================================================
	samples = []
	sequence_files = []
	mysplit = args.split_file
	# myfile = os.path.join(args.input, args.prefix + "." + mysplit + ".fasta")
	myfile = args.input
	samples.append(mysplit)
	sequence_files.append(myfile)

	# open_file = open(args.file_list, "r")
	# for line in open_file.readlines():
	#    line = line.strip()
	#    if not len(line):
	#       continue
	#    mysplit = line
	#    myfile = os.path.join(args.input, mysplit, args.prefix + "." + mysplit + ".fasta")
	#    samples.append(mysplit)
	#    sequence_files.append(myfile)
	# foreach file
	# open_file.close()

	## PSORTb: protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes
	annotation_dir = args.output

	for protein in sequence_files:
		# Using Gram+
		protein_base = os.path.basename(protein).split(os.extsep)[-2]
		out_file = os.path.join(annotation_dir, '%s.psortb.gram_positive.out.txt' % protein_base)
		stderr_log = os.path.join(annotation_dir, '%s.psortb.gram_positive.err' % protein_base)
		os.system("psort -p " + protein + " > " + out_file + " 2> " + stderr_log)

		# Gram-
		out_file = os.path.join(annotation_dir, '%s.psortb.gram_negative.out.txt' % protein_base)
		stderr_log = os.path.join(annotation_dir, '%s.psortb.gram_negative.err' % protein_base)
		os.system("psort -n " + protein + " > " + out_file + " 2> " + stderr_log)

		# Archaea
		out_file = os.path.join(annotation_dir, '%s.psortb.archaea.out.txt' % protein_base)
		stderr_log = os.path.join(annotation_dir, '%s.psortb.archaea.err' % protein_base)
		os.system("psort -a " + protein + " > " + out_file + " 2> " + stderr_log)
	# foreach sequence file



if __name__ == "__main__":
	main()

#!/usr/bin/env python

"""
interproscan_annotator.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run interproscan for protein annotations

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

import os
import sys 
import argparse
import csv 
import re

try:
	from metawibele import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

description = """
A workflow to annotate domain signatures
"""

def parse_cli_arguments():
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
	# sequence_files = glob(os.path.join(args.input, "*%s" % args.file_extension))
	# samples = [os.path.basename(s).split(os.extsep)[0] for s in sequence_files]
	samples = []
	sequence_files = []
	mysplit = args.split_file
	myfile = args.input
	samples.append(mysplit)
	sequence_files.append(myfile)

	## InterProScan 5: genome-scale protein function classification
	annotation_dir = args.output

	for protein in sequence_files:
		protein_base = os.path.basename(protein).split(os.extsep)[-2]
		out_file = os.path.join(annotation_dir, '%s.interproscan.txt' % protein_base)
		stderr_log = os.path.join(annotation_dir, '%s.interproscan.err' % protein_base)
		#commands = config.interproscan_cmmd + " -appl " + config.interproscan_appl + " -cpu " + str(args.threads) + " -i " + protein + " -f tsv -dp -t p -o " + out_file + " >" + stderr_log + " 2>&1"
		commands = config.interproscan_cmmd + " -appl " + config.interproscan_appl + " -i " + protein + " -f tsv -dp -t p -o " + out_file + " >" + stderr_log + " 2>&1"
		print(commands)
		os.system(commands)
	# foreach sequence file


if __name__ == "__main__":
	main()

#!/usr/bin/env python

"""
MetaWIBELE: download_config_file module
Download template configuration files and vignette files for prioritization used by MetaWIBELE

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

# Try to load one of the MetaWIBELE src modules to check the installation
try:
	from metawibele import config, check
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")


# Check the python version
check.python_version()

import argparse

config_types={
    "global" : "metawibele.cfg",
    "local" : "characterization.cfg,prioritization.cfg,MSPminer_setting.cfg",
	"vignette": "vignettes_proteins.txt.gz"
}

def download_config(config_type, config_types):
	"""
	Download the selected configuration files
	"""
	if config_type in config_types:
		files = config_types[config_type].split(",")
		print("Downloading template configuration files...")
		for i in files:
			file_location = os.path.join(config.config_directory, i)
			if os.path.isfile(file_location):
				os.system("cp " + file_location + " " + os.getcwd())
			else:
				file_location = os.path.join(config.misc_directory, i)
				if os.path.isfile(file_location):
					os.system("cp " + file_location + " " + os.getcwd())
				else:
					print("The configuration file selected does not exist! " + file_location)

		print("\nFinished downloading configuration files!\n")
	else:
		sys.exit("ERROR: Please select an available type of configuration file.")

def parse_arguments(args):
	"""
	Parse the arguments from the user
	"""
	parser = argparse.ArgumentParser(
		description = "Download MetaWIBELE configuration files\n",
		formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument(
		"--config-type",
		choices=["global", "local", "vignette"],
		required=True)

	return parser.parse_args()

def main():
	# Parse arguments from the command line
	args = parse_arguments(sys.argv)

	download_config(args.config_type, config_types)

if __name__ == '__main__':
	main()

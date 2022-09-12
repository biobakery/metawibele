#!/usr/bin/env python

"""
MetaWIBELE: check_install module
Check the installation of MetaWIBELE package and dependent tools

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
import datetime
import time
import re
import subprocess
import argparse
from pathlib import Path

def check_metawibele ():
	"""
	Check metawibele package
	"""

	print("###### Checking anadama2 package install ######")
	try:
		import anadama2
		print("Checking anadama2 install...........OK")
	except ImportError:
		print("CRITICAL ERROR: Unable to find anadama2 package. " +
		      "Please check your install.")

	print("\n###### Checking metawibele package install ######")
	try:
		import metawibele
		print("Checking metawibele install...........OK")
		from metawibele import config
		if config.uniref_database_dir.lower() == "none" or config.uniref_database_dir == "":
			pass
			#print("\n  WARNING!! MetaWIBELE does't find valid uniref databse and will use the demo databse by default.\n" +
			#	"\tPlease provide the correct location of the required uniref database by any one of the following options:\n" +
			#	"\t1) set the location with the environment variable $UNIREF_LOCATION\n" +
			#	"\t2) include the database (named as \"uniref_database\") in the current working directory\n" +
			#	"\t3) set the location in the global config file (metawibele.cfg) which is in the current working directory"
			#      )

	except ImportError:
		print("CRITICAL ERROR: Unable to find the MetaWIBELE python package. " +
		         "Please check your install.")


def run_cmd (cmd):
	out = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
	mystdout, mystderr = out.communicate()
	if mystdout:
		mystdout = mystdout.decode()
	if mystderr:
		mystderr = mystderr.decode()
	mystderr = None
	return mystdout, mystderr


def check_required_tools ():
	"""
	Check required dependent tools
	"""

	print("\n###### Checking CD-hit install ######")
	try:
		mystout, mysterr = run_cmd (['cd-hit', '-h'])
		if mysterr:
			print("Checking CD-hit install...........Error!")
			print(mysterr)
		else:
			print("Checking CD-hit install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking CD-hit install...........Error!")
		os.system("cd-hit -h")

	print("\n###### Checking Diamond install ######")
	try:
		mystout, mysterr = run_cmd (['diamond', '--version'])
		if mysterr:
			print("Checking Diamond install...........Error!")
			print(mysterr)
		else:
			print("Checking Diamond install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking Diamond install...........Error!")
		os.system("diamond --version")

	print("\n###### Checking MSPminer install ######")
	try:
		mystout, mysterr = run_cmd (['mspminer'])
		if mysterr:
			print("Checking MSPminer install...........Error!")
			print(mysterr)
		else:
			print("Checking MSPminer install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking MSPminer install...........Error!")
		os.system("mspminer")

	print("\n###### Checking Maaslin2 install ######")
	try:
		from metawibele import config
		maaslin_cmd = config.maaslin2_cmmd
		if re.search("^\~", maaslin_cmd):
			home = str(Path.home())
			maaslin_cmd = re.sub("^~", home, maaslin_cmd)
	except ImportError:
		maaslin_cmd = "Maaslin2.R"
	try:
		mystout, mysterr = run_cmd ([maaslin_cmd, '-h'])
		if mysterr:
			print("Checking MaAsLin2 install...........Error! Please check the valid location for " + maaslin_cmd)
			print(mysterr)
		else:
			print("Checking MaAsLin2 install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking MaAsLin2 install...........Error! Please check the valid location for " + maaslin_cmd)
		os.system(maaslin_cmd + " -h")

	print("\n###### Checking Interproscan install ######")
	try:
		from metawibele import config
		interproscan_cmd = config.interproscan_cmmd
		if re.search("^\~", interproscan_cmd):
			home = str(Path.home())
			interproscan_cmd = re.sub("^~", home, interproscan_cmd)
	except ImportError:
		maaslin_cmd = "interproscan.sh"
	try:
		mystout, mysterr = run_cmd ([interproscan_cmd, '-version'])
		if mysterr:
			print("Checking Interproscan install...........Error! Please check valid location for " + interproscan_cmd + "\n" +
			      "\tPlease install Interproscan in a location in your $PATH or provide the location to MetaWIBELE global config file.\n" +
			      "\tIf you'd like to skip annotate proteins based on domain/motif, use the option '--bypass-domain-motif' or '--bypass-interproscan' when running MetaWIBELE."
			      )
			print(mysterr)
		else:
			print("Checking Interproscan install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking Interproscan install...........Error! Please check valid location for " + interproscan_cmd + "\n" +
		      "\tPlease install Interproscan in a location in your $PATH or provide the location to MetaWIBELE global config file.\n" +
		      "\tIf you'd like to skip annotate proteins based on domain/motif, use the option '--bypass-domain-motif' or '--bypass-interproscan' when running MetaWIBELE."
		      )
		os.system(interproscan_cmd + " -version")

	print("\n###### Checking Signalp install ######")
	try:
		mystout, mysterr = run_cmd (['signalp', '-V'])
		if mysterr:
			print("Checking Signalp install...........Warning!\n" +
			      "\tSignalp is not installed in a location in your $PATH. Please install it and provide the install location to interproscan.properties to active signalp analysis when running interproscan.\n" +
				  "\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			      )
			print(mysterr)
		else:
			print("Checking Signalp install...........Inconclusive\n" +
			      "\tPlease make sure the install location of Signalp has been provided to interproscan.properties for activing signalp analysis when running interproscan.\n" +
			      "\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			      )
			if mystout:
				print(mystout)
	except:
		print("Checking Signalp install...........Warning!\n" +
			"\tSignalp is not installed in a location in your $PATH. Please install it and provide the install location to interproscan.properties to active signalp analysis when running interproscan.\n" +
			"\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			)
		os.system("signalp -V") 

	print("\n###### Checking TMHMM install ######")
	try:
		mystout, mysterr = run_cmd (['which', 'tmhmm'])
		if mysterr:
			print("Checking TMHMM install...........Warning!\n" +
			      "\tTMHMM is not installed in a location in your $PATH. Please install it and provide the install location to interproscan.properties to active TMHMM analysis when running interproscan.\n" +
				  "\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			      )
			print(mysterr)
		else:
			print("Checking TMHMM install...........Inconclusive\n" +
			      "\tPlease make sure the install location of TMHMM has been provided to interproscan.properties for activing TMHMM analysis when running interproscan.\n" +
			      "\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			      )
			if mystout:
				print(mystout)
	except:
		print("Checking TMHMM install...........Warning!\n" +
			"\tTMHMM is not installed in a location in your $PATH. Please install it and provide the install location to interproscan.properties to active TMHMM analysis when running interproscan.\n" +
			"\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			)
		os.system("which tmhmm")

	print("\n###### Checking Phobius install ######")
	try:
		mystout, mysterr = run_cmd (['phobius.pl', '-h'])
		if mysterr:
			print("Checking Phobius install...........Warning!\n" +
			      "\tPhobius is not installed in a location in your $PATH. Please install it and provide the install location to interproscan.properties to active Phobius analysis when running interproscan.\n" +
				  "\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			      )
			print(mysterr)
		else:
			print("Checking Phobius install...........Inconclusive\n" +
			      "\tPlease make sure the install location of Phobius has been provided to interproscan.properties for activing Phobius analysis when running interproscan.\n" +
			      "\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
				)
			if mystout:
				print(mystout)
	except:
		print("Checking Phobius install...........Warning!\n" +
			"\tPhobius is not installed in a location in your $PATH. Please install it and provide the install location to interproscan.properties to active Phobius analysis when running interproscan.\n" +
			"\tSee more details from InterProScan document (https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)."
			)
		os.system("phobius.pl -h")

	print("\n###### Checking PSORTb install ######")
	try:
		mystout, mysterr = run_cmd (['psort', '--version'])
		if mysterr:
			print("Checking PSORTb install...........Error!\n" +
			      "\tPlease install PSORTb in a location in your $PATH.\n" +
			      "\tIf you'd like to skip annotate proteins based on domain/motif by PSORTb, use the option '--bypass-psortb' when running MetaWIBELE."
			      )
			print(mysterr)
		else:
			print("Checking PSORTb install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking PSORTb install...........Error!\n" +
		      "\tPlease install PSORTb in a location in your $PATH.\n" +
		      "\tIf you'd like to skip annotate proteins based on domain/motif by PSORTb, use the option '--bypass-psortb' when running MetaWIBELE."
		      )
		os.system("psort --version")


def check_optional_tools ():
	"""
	Check optional dependent tools for using MetaWIBELE utilities to prepare MetaWIBELE inputs
	"""
	print("\n###### Checking MEGAHIT install ######")
	try:
		mystout, mysterr = run_cmd (['megahit', '-v'])
		if mysterr:
			print("Checking MEGAHIT install...........Error!\n" +
			      "\tPlease install MEGAHIT in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			      )
			print(mysterr)
		else:
			print("Checking MEGAHIT install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking MEGAHIT install...........Error!\n" +
		      "\tPlease install MEGAHIT in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			)
		os.system("megahit -v")

	print("\n###### Checking Prokka install ######")
	try:
		mystout, mysterr = run_cmd (['prokka', '-v'])
		if mysterr:
			print("Checking Prokka install...........Error!\n" +
			      "\tPlease install Prokka in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			      )
			print(mysterr)
		else:
			print("Checking Prokka install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking Prokka install...........Error!\n" +
		      "\tPlease install Prokka in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			)
		os.system("prokka -v")

	print("\n###### Checking Prodigal install ######")
	try:
		mystout, mysterr = run_cmd (['prodigal', '-v'])
		if mysterr:
			print("Checking Prodigal install...........Error!\n" +
			      "\tPlease install Prodigal in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			      )
			print(mysterr)
		else:
			print("Checking Prodigal install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking Prodigal install...........Error!\n" +
		      "\tPlease install Prodigal in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			)
		os.system("prodigal -v")

	print("\n###### Checking USEARCH install ######")
	try:
		mystout, mysterr = run_cmd (['usearch', '--help'])
		if mysterr:
			print("Checking USEARCH install...........Error!\n" +
			      "\tPlease install USEARCH in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			      )
			print(mysterr)
		else:
			print("Checking USEARCH install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking USEARCH install...........Error!\n" +
		      "\tPlease install USEARCH in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			)
		os.system("usearch --help")

	print("\n###### Checking Bowtie2 install ######")
	try:
		mystout, mysterr = run_cmd (['bowtie2', '--version'])
		if mysterr:
			print("Checking Bowtie2 install...........Error!\n" +
			      "\tPlease install Bowtie2 in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			      )
			print(mysterr)
		else:
			print("Checking Bowtie2 install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking Bowtie2 install...........Error!\n" +
		      "\tPlease install Bowtie2 in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			)
		os.system("bowtie2 --version")

	print("\n###### Checking SAMtools install ######")
	try:
		mystout, mysterr = run_cmd (['samtools', '--version'])
		if mysterr:
			print("Checking SAMtools install...........Error!\n" +
			      "\tPlease install SAMtools in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			      )
			print(mysterr)
		else:
			print("Checking SAMtools install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking SAMtools install...........Error!\n" +
		      "\tPlease install SAMtools in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			)
		os.system("samtools --version")

	print("\n###### Checking featureCounts install ######")
	try:
		mystout, mysterr = run_cmd (['featureCounts', '-v'])
		if mysterr:
			print("Checking featureCounts install...........Error!\n" +
			      "\tPlease install featureCounts in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			      )
			print(mysterr)
		else:
			print("Checking featureCounts install...........OK")
			if mystout:
				print(mystout)
	except:
		print("Checking featureCounts install...........Error!\n" +
		      "\tPlease install featureCounts in in a location in your $PATH if you will use MetaWIBELE utilities to prepare MetaWIBELE inputs."
			)
		os.system("featureCounts -v")


def parse_arguments(args):
	"""
	Parse the arguments from the user
	"""
	parser = argparse.ArgumentParser(
		description = "Check the installation of MetaWIBELE package and dependent tools\n")
	parser.add_argument(
		"--types",
		help = "provide which type of tools for checking",
		choices = ["metawibele", "required", "optional", "all"],
		default = "all")

	return parser.parse_args()


def main():
	# Parse arguments from the command line
	args = parse_arguments(sys.argv)

	if args.types == "metawibele":
		print("#----------------- Checking the install of metawibele packages -----------------#")
		check_metawibele()

	if args.types == "required":
		print("#----------------- Checking the install of required dependent tools -----------------#")
		check_required_tools()

	if args.types == "optional":
		print("#----------------- Checking the install of optional dependent tools -----------------#")
		check_optional_tools()

	if args.types == "all":
		print("#----------------- Checking the install of metawibele packages -----------------#")
		check_metawibele()
		print("\n\n#----------------- Checking the install of required dependent tools -----------------#")
		check_required_tools()
		print("\n\n#----------------- Checking the install of optional dependent tools -----------------#")
		check_optional_tools()

if __name__ == '__main__':
	main()

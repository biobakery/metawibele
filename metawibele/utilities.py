#!/usr/bin/env python
"""
MetaWIBELE: utilities module
Utilities relating to third party software, file permissions, and file formats

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
import os.path
import re
import argparse
import subprocess
import csv
import gzip
import bz2
import time
import math

#from anadama2.tracked import TrackedDirectory

# constants
c_strat_delim = "|"
c_taxon_delim = "."
c_name_delim = "__"
c_multiname_delim = ";"
c_str_unknown = "Unknown"
c_unclassified = "Unclassified"
c_many_bytes = 1e8
c_zip_multiplier = 10

PROTEIN_FAMILY_ID = "familyID"
PROTEIN_ID = "seqID"
BIOM_FILE_EXTENSION = ".biom"  # the extension used for biom files


# ==============================================================
# utilities used for metadata, cluster, GO info collection
# ==============================================================
def sample_info (sampleinfo, study):
	""" Collect metadata info and return sample info"""

	samples = {}
	title = {}
	if not os.path.isfile(sampleinfo):
		print ("File not exist!\t" + sampleinfo)
		return
	open_file = open(sampleinfo, "r")
	line = open_file.readline()
	line = line.strip()
	info = line.split("\t")
	for item in info:
		title[item] = info.index(item)
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if "External_ID" in title:
			sample = info[title["External_ID"]]
		if "ID" in title:
			sample = info[title["ID"]]
		if "SID" in title:
			sample = info[title["SID"]]
		disease = info[title["diagnosis"]]
		if study == "HMP2":
			disease = info[title["diagnosis2"]]
		if re.search("^[\d]+$", sample):
			sample = study + "_" + sample
		sample = re.sub("^" + disease + "_", "", sample)
		samples[sample] = disease
	# foreac sample
	open_file.close()
	return samples


# sample_info


def collect_partial_info (infile):
	"""
	collect partial info
	Input: gene details file
	Output: partial gene info file
	"""

	partial = {}
	titles = {}
	open_file = open(infile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^sample", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		gene = info[titles["gene"]]
		mypartial = info[titles["partial"]]
		feature = info[titles["feature"]]
		start = info[titles["start"]]
		stop = info[titles["end"]]
		mylen = abs(int(stop) - int(start) + 1)
		if feature != "CDS":
			continue
		flag = "NA"
		if mypartial == "00":  # complete genes
			flag = "complete"
		if mypartial == "10":  # lose start codon
			flag = "no_start_codon"
		if mypartial == "01":  # lose stop codon
			flag = "no_stop_codon"
		if mypartial == "11":  # no stop and start
			flag = "no_start_stop"
		partial[gene] = flag
	# foreach line
	open_file.close()
	return partial


# collect_partial_info


def collect_partial_info_cluster (clust_file):
	"""
	Collect partial info for clusters
	Input: partial gene info within protein families, PRISM_clustering.Linclust.complete_ORF.detail.tsv
	Output: partial info for protein families
	"""

	partial = {}
	titles = {}
	open_file = open(clust_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^cluster\t", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myclust = info[titles["cluster"]]  # Cluster_1
		mytype = info[titles["type"]]
		partial[myclust] = mytype
	# foreach line
	return partial


# collect_partial_info_cluster


def collect_partial_all_info_cluster (clust_file):
	"""
	Collect partial info for clusters
	Input: partial gene info within protein families, PRISM_clustering.Linclust.complete_ORF.detail.tsv
	Output: partial info for protein families
	"""

	partial = {}
	titles = {}
	open_file = open(clust_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if re.search("^cluster\t", line):
			for item in info:
				titles[item] = info.index(item)
			continue
		myclust = info[titles["cluster"]]  # Cluster_1
		mytype = info[titles["type"]]
		mycat = info[titles["category"]]
		partial[myclust] = mytype + "\t" + mycat
	# foreach line
	return partial


# collect_partial_all_info_cluster


def collect_protein_cluster_info (clust_file):
	"""
	Collect protein family info
	Return: {centroid id: {gene id : cluster id}}
	"""

	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	myclust_id = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			mym = re.search("cluster=([\d]+)", line)
			myclust_id = "Cluster_" + mym.group(1)
			if not myclust in cluster:
				cluster[myclust] = {}
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myclust][myid] = myclust_id
	# foreach line
	open_file.close()
	return cluster

# function collect_protein_cluster_info


def collect_gene_cluster_info (clust_file):
	"""
	Collect gene catalog info
	Return: {gene id : cluster id}
	"""

	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search(">([^;]+)", line)
			myclust = mym.group(1)
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myclust
	# foreach line
	open_file.close()
	return cluster

# function collect_gene_cluster_info


def collect_GO_info(gofile):  # go.obo
	"""
	Collect GO term info
	"""
	gos = {}
	open_file = open(gofile, "r")
	mygo = "NA"
	myname = "NA"
	mytype = "NA"
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("id: ", line):
			mym = re.search("id: ([\S]+)", line)
			mygo = mym.group(1)
		if re.search("name: ", line):
			mym = re.search("name: ([\S\s]+)", line)
			myname = mym.group(1)
		if re.search("^namespace: ", line):
			mym = re.search("namespace: ([\S]+)", line)
			mytype = mym.group(1)
			if mytype == "biological_process":
				mytype = "GO(BP)"
			if mytype == "molecular_function":
				mytype = "GO(MF)"
			if mytype == "cellular_component":
				mytype = "GO(CC)"
			gos[mygo] = myname + "[" + mygo + "]" + "\t" + mytype
	# foreach line
	open_file.close()
	return gos


# collect_GO_info


# ==============================================================
# utilities used for basis statistics info
# ==============================================================
def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data) / n  # in Python 2 use sum(data)/float(n)


def _ss(data):
	"""Return sum of square deviations of sequence data."""
	c = mean(data)
	ss = sum((x - c) ** 2 for x in data)
	return ss


def stddev(data, ddof=0):
	"""Calculates the population standard deviation
	by default; specify ddof=1 to compute the sample
	standard deviation."""
	n = len(data)
	if n < 2:
		raise ValueError('variance requires at least two data points')
	ss = _ss(data)
	pvar = ss / (n - ddof)
	return pvar ** 0.5


def remove_duplicate(duplicate):
	"""
	remove duplicates in list
	"""
	final_list = []
	for num in duplicate:
		if num not in final_list:
			final_list.append(num)
	return final_list


# remove_duplicate


# ==============================================================
# utilities used for run tasks
# =============================================================
def run_task(command, **keywords):
	""" Run the task command, formatting command with keywords. The command stdout
		and stderr are written to the workflow log.

	Args:
		command (string): A string to execute on the command line. It can be
			formatted the same as a task command.

	Returns:
		(int): Return code from command.
	"""

	from anadama2.helpers import format_command
	from anadama2.helpers import sh

	# format the command to include the items for this task
	command = format_command(command, **keywords)

	# run the command
	return_code = sh(command)()

	return return_code


# ==============================================================
# utilities used for I/O files
# =============================================================
def find_files(folder, extension=None, exit_if_not_found=None):
	""" Return the files in the given folder with the extension if provided

	Args:
		folder (string): A path to a folder
		extension (string): The file extension to search for (optional)
		exit_if_not_found (bool): Indicator to check if files exist (optional)

	Requires:
		None

	Returns:
		list: A list of files in the folder

	Example:
		files = find_files("examples","fastq")
	"""

	# get all of the files in the folder
	files = [os.path.join(folder, file) for file in os.listdir(folder)]
	files = list(filter(lambda file: os.path.isfile(file), files))

	# filter to only files with extension
	if extension:
		files = list(filter(lambda file: file.endswith(extension), files))

	if exit_if_not_found:
		if not files:
			message = "ERROR: No files were found in the folder " + folder
			if extension:
				message += " with extension " + extension
			sys.exit(message + " .\n")

	return files


def name_files(names, folder, subfolder=None, tag=None, extension=None, create_folder=None):
	""" Return a list of file names based on the names and folders provided

	Args:
		names (list or string): A list of basenames or files.
		folder (string): The path to the folder.
		subfolder (string): The subfolder to use with the files (optional).
		tag (string): The tag to add to the file basenames (optional).
		extension (string): The extension to use for the files (optional).
		create_folder (bool): Create the folder and subfolder if they do not exist (optional).

	Requires:
		None

	Returns:
		list: A list of file names (or string if input is string).

	Example:
		files = name_files(["file1","file2"], "output")
	"""

	# if names is a list, convert to string
	was_string = False
	if isinstance(names, basestring):
		was_string = True
		names = [names]

	# get the basenames from the files
	names = [os.path.basename(name) for name in names]

	# use the full path to the folder
	folder = os.path.abspath(folder)

	# get the name of the full folder plus subfolder if provided
	if subfolder:
		folder = os.path.join(folder, subfolder)

	# add the extension if provided, and replace existing
	if extension:
		names = [os.path.splitext(name)[0] + "." + extension for name in names]

	# add the tag to the names, if provided
	if tag:
		names = [os.path.splitext(name)[0] + "_" + tag + os.path.splitext(name)[1] for name in names]

	files = [os.path.join(folder, name) for name in names]

	if create_folder:
		create_folders(os.path.dirname(files[0]))

	# if the input was originally a string, convert from list
	if was_string:
		files = files[0]

	return files


def create_folders(folder):
	""" Create folder if it does not exist

	Args:
		folder (string): The full path to the folder.

	Requires:
		None

	Returns:
		None

	Example:
		create_folders("new_folder")
	"""

	try:
		if not os.path.exists(folder):
			os.makedirs(folder)
	except EnvironmentError:
		print("Warning: Unable to create folder: " + folder)


def get_files(folder, extension):
	""" Return paths to all files in a folder with a given extension """

	for file in os.listdir(folder):
		file = os.path.join(folder, file)
		if os.path.isfile(file) and file.endswith(extension):
			yield file


def find_exe_in_path(exe):
	"""
	Check that an executable exists in $PATH
	"""

	paths = os.environ["PATH"].split(os.pathsep)
	for path in paths:
		fullexe = os.path.join(path, exe)
		if os.path.exists(fullexe):
			if os.access(fullexe, os.X_OK):
				return True
	return False


def biom_to_tsv(biom_file, new_tsv_file, taxonomy=None):
	"""
	Convert from a biom to tsv file
	"""

	cmd = ["biom", "convert", "-i", biom_file, "-o", new_tsv_file, "--to-tsv"]

	# check if taxonomy is set (can be set to zero)
	if taxonomy != None:
		cmd += ["--header-key", "taxonomy"]

	try:
		if os.path.isfile(new_tsv_file):
			os.unlink(new_tsv_file)
		p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
	except (EnvironmentError, subprocess.CalledProcessError):
		command = " ".join(cmd)
		sys.exit("Unable to convert biom file to tsv" + "\n" + command)


def tsv_to_biom(tsv_file, biom_file):
	"""
	Convert from a biom to tsv file
	"""

	cmd = ["biom", "convert", "-i", tsv_file, "-o", biom_file, "--table-type", "Gene table", "--to-hdf5"]

	try:
		if os.path.isfile(biom_file):
			os.unlink(biom_file)
		p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
	except (EnvironmentError, subprocess.CalledProcessError):
		command = " ".join(cmd)
		sys.exit("Unable to convert tsv file to biom" + "\n" + command)


def process_gene_table_with_header(gene_table, allow_for_missing_header=None):
	"""
	Process through the header portion of the gene table file
	"""

	# try to open the file
	try:
		lines = gzip_bzip2_biom_open_readlines(gene_table)
	except EnvironmentError:
		sys.exit("Unable to read file: " + gene_table)

	# find the headers
	header = ""
	first_data_line = ""
	for line in lines:
		if line[0] == GENE_TABLE_COMMENT_LINE:
			header = line
		else:
			first_data_line = line
			break

	if not header and not allow_for_missing_header:
		sys.exit("File does not have a required header: " + gene_table +
		         " . Please add a header which includes the indicator: " +
		         GENE_TABLE_COMMENT_LINE)

	# provide the header, if one was found
	if header:
		yield header

	# provide the first data line
	yield first_data_line

	# now provide the remaining lines
	for line in lines:
		yield line


def write_tsv(path, rows):
	""" Write the output in tsv (possibly compressed) format to a file or stdout """
	fh = try_zip_open(path, write=True) if path is not None else sys.stdout
	writer = csv.writer(fh, delimiter="\t", lineterminator="\n")

	for row in rows:
		writer.writerow(row)


def write_biom(path, rows):
	""" Write the file in biom format """

	try:
		import biom
	except ImportError:
		sys.exit("Could not find the biom software." +
		         " This software is required since the input file is a biom file.")

	try:
		import numpy
	except ImportError:
		sys.exit("Could not find the numpy software." +
		         " This software is required since the input file is a biom file.")

	try:
		import h5py
	except ImportError:
		sys.exit("Could not find the h5py software." +
		         " This software is required since the input file is a biom file.")

	# reformat the rows into a biom table
	samples = next(rows)[1:]
	ids = []
	data = []
	for row in rows:
		ids.append(row[0])
		data.append(row[1:])

	table = biom.Table(numpy.array(data), ids, samples)

	# write a h5py biom table
	with h5py.File(path, 'w') as file_handle:
		table.to_hdf5(file_handle, "metawibele utility script")


class Ticker():
	def __init__(self, iterable, step=100, pad="  "):
		self.count = 0
		self.total = len(iterable)
		self.step = 100
		self.pad = pad

	def tick(self):
		self.count += 1
		if self.count % self.step == 0:
			self.report()

	def report(self):
		frac = self.count / float(self.total)
		# print(self.pad + "{:.1f}%".format(100 * frac), file=sys.stderr, end="\r")


# ---------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------

def size_warn(path):
	m = 1 if ".gz" not in path else c_zip_multiplier
	if m * os.path.getsize(path) > c_many_bytes:
		# print("  This is a large file, one moment please...", file=sys.stderr)
		print("  This is a large file, one moment please...")


def try_zip_open(path, write=None):
	"""
	open an uncompressed or gzipped file; fail gracefully
	"""
	fh = None

	# set the open mode
	if write:
		open_mode = "w"
	elif path.endswith(".bz2"):
		open_mode = "r"
	else:
		open_mode = "rt"

	try:
		if path.endswith(".gz"):
			fh = gzip.open(path, open_mode)
		elif path.endswith(".bz2"):
			fh = bz2.BZ2File(path, open_mode)
		else:
			fh = open(path, open_mode)
	except EnvironmentError:
		sys.exit("Problem opening file: " + path)
	return fh


def read_biom_table(path):
	"""
	return the lines in the biom file
	"""

	try:
		import biom
	except ImportError:
		sys.exit("Could not find the biom software." +
		         " This software is required since the input file is a biom file.")

	try:
		tsv_table = biom.load_table(path).to_tsv().split("\n")
	except (EnvironmentError, TypeError):
		sys.exit("ERROR: Unable to read biom input file.")

	return tsv_table


def gzip_bzip2_biom_open_readlines(path):
	"""
	return the lines in the opened file for tab delimited text, gzip, bzip2 and biom files
	"""

	# if the file is biom, convert to text and return lines
	if path.endswith(BIOM_FILE_EXTENSION):
		for line in read_biom_table(path):
			yield line
	else:
		with try_zip_open(path) as file_handle:
			for line in file_handle:
				if path.endswith(".bz2"):
					# convert the line to text from binary
					yield line.decode('utf-8').rstrip()
				else:
					yield line.rstrip()


def fsplit(feature):
	items = feature.split(c_strat_delim)
	stratum = None if len(items) == 1 else items[1]
	items = items[0].split(c_name_delim)
	name = None if len(items) == 1 else items[1]
	feature = items[0]
	return feature, name, stratum


def fjoin(feature, name=None, stratum=None):
	if name is not None:
		feature = c_name_delim.join([feature, name])
	if stratum is not None:
		feature = c_strat_delim.join([feature, stratum])
	return feature


def fsort(features):
	# force 1|A to come before 11
	features = sorted(features, key=lambda f: f.split(c_strat_delim))
	# force special features to the top (defined above)
	default = 1 + max(c_topsort.values())
	features = sorted(features, key=lambda f: c_topsort.get(fsplit(f)[0], default))
	return features


if __name__ == '__main__':
	pass
#!/usr/bin/env python

"""
MetaWIBELE: download_databases module
Download dependent databases used by MetaWIBELE

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

try:
	from urllib.request import urlretrieve
except ImportError:
	from urllib import urlretrieve
import tarfile
import gzip

# Try to load one of the MetaWIBELE src modules to check the installation
try:
	from metawibele import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
	         " Please check your install.")

import argparse

# the locations of the current databases to download
current_downloads = {
	"uniref":
		{
			"uniref90_diamond": "https://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90.fasta.dmnd.tar.gz",
			"uniref90_fasta": "https://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90.fasta.tar.gz",
			"uniref90_annotation": "https://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90_annotations.tar.gz"
		}
}


def byte_to_gigabyte(byte):
	"""
	Convert byte value to gigabyte
	"""

	return byte / (1024.0 ** 3)


def byte_to_megabyte(byte):
	"""
	Convert byte value to megabyte
	"""
	
	return byte / (1024.0**2)


def byte_to_kilobyte(byte):
	"""
	Convert byte value to kilobyte
	"""
	
	return byte / 1024.0


class ReportHook():
	def __init__(self):
		self.start_time = time.time()

	def report(self, blocknum, block_size, total_size):
		"""
		Print download progress message
		"""

		if blocknum == 0:
			self.start_time = time.time()
			if total_size > 0:
				print("Downloading file of size: " + "{:.2f}".format(byte_to_gigabyte(total_size)) + " GB\n")
		else:
			total_downloaded = blocknum * block_size
			status = "{:3.2f} GB ".format(byte_to_gigabyte(total_downloaded))

			if total_size > 0:
				percent_downloaded = total_downloaded * 100.0 / total_size
				# use carriage return plus sys.stdout to overwrite stdout
				download_rate = total_downloaded / (time.time() - self.start_time)
				estimated_time = (total_size - total_downloaded) / download_rate
				estimated_minutes = int(estimated_time / 60.0)
				estimated_seconds = estimated_time - estimated_minutes * 60.0
				status += "{:3.2f}".format(percent_downloaded) + " %  " + \
				          "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
				          "{:2.0f}".format(estimated_minutes) + " min " + \
				          "{:2.0f}".format(estimated_seconds) + " sec "
			status += "        \r"
			sys.stdout.write(status)


def download_file(url, filename, folder):
	"""
	Download the file at the url
	"""

	print("Download URL: " + url)

	try:
		#url_handle = urlretrieve(url, filename, reporthook=ReportHook().report)
		os.system("wget " + url + " --no-check-certificate ")
		print("\nExtracting: " + filename)
		os.system("tar zxf " + filename)
		#tarfile_handle = tarfile.open(filename)
		#tarfile_handle.extractall(path=folder)
	except (EnvironmentError, tarfile.ReadError):
		sys.exit("CRITICAL ERROR: Unable to download and extract from URL: " + url)


def download_database(database, db_type, location):
	"""
	Download and decompress the selected database
	"""

	if database in current_downloads:
		if db_type in current_downloads[database]:
			if not os.path.isdir(location):
				try:
					print("Creating directory to install database: " + location)
					os.mkdir(location)
				except EnvironmentError:
					sys.exit("CRITICAL ERROR: Unable to create directory: " + location)

			# download the database
			downloaded_file = os.path.join(location, current_downloads[database][db_type].split('/')[-1])
			download_file(current_downloads[database][db_type], downloaded_file, location)

			# remove the download
			try:
				os.unlink(downloaded_file)
			except EnvironmentError:
				print("Unable to remove file: " + downloaded_file)

			print("\nDatabase installed: " + location + "\n")
		else:
			sys.exit("ERROR: Please select an available database type.")
	else:
		sys.exit("ERROR: Please select an available database.")


def parse_arguments(args):
	"""
	Parse the arguments from the user
	"""
	parser = argparse.ArgumentParser(
		description="Download MetaWIBELE dependent databases\n",
		formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument(
		"--database",
		help="provide which database do you want to download",
		choices=["uniref"],
		default="uniref")
	parser.add_argument(
		"--build",
		help="provide which type of data do you want to download from the database",
		choices=["uniref90_diamond", "uniref90_fasta", "uniref90_annotation"],
		required=True)
	parser.add_argument(
		"--install-location",
		help="provide the location for installing the database",
		required=True)

	return parser.parse_args()


def main():
	# Parse arguments from the command line
	args = parse_arguments(sys.argv)

	download_database(args.database, args.build, args.install_location)


if __name__ == '__main__':
	main()

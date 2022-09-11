"""
MetaWIBELE setup

To run: python setup.py install
"""

import os
import sys

# check for either of the required versions
# required python versions (3.6+)
required_python_version_major = [3]
required_python_version_minor = [6]
pass_check = False
try:
	for major, minor in zip(required_python_version_major, required_python_version_minor):
		if (sys.version_info[0] == major and sys.version_info[1] >= minor):
			pass_check = True
except (AttributeError, IndexError):
	sys.exit("CRITICAL ERROR: The python version found (version 1) " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")


try:
	import setuptools
except ImportError:
	sys.exit("Please install setuptools.")

# check setuptools version
required_setuptools_version_major = 1
try:
	setuptools_version = setuptools.__version__
	setuptools_version_major = int(setuptools_version.split(".")[0])
	if setuptools_version_major < required_setuptools_version_major:
		sys.exit("CRITICAL ERROR: The setuptools version found (version " +
		         setuptools_version + ") does not match the version required " +
		         "(version " + str(required_setuptools_version_major) + "+)." +
					" Please upgrade your setuptools version.")
except (ValueError, IndexError, NameError):
	sys.exit("CRITICAL ERROR: Unable to call setuptools version. Please upgrade setuptools.")

from setuptools.command.install import install as _install

import distutils

# try to import urllib.request.urlretrieve for python3
try:
	from urllib.request import urlretrieve
except ImportError:
	from urllib import urlretrieve

from glob import glob
import re
import tarfile
import subprocess
import shutil
import zipfile
import tempfile
import re
import time


VERSION = "0.4.6"
AUTHOR = "MetaWIBELE Development Team"
MAINTAINER = "Yancong Zhang"
MAINTAINER_EMAIL = "zhangyc201211@gmail.com"


def byte_to_megabyte(byte):
	"""
	Convert byte value to megabyte
	"""

	return byte / (1024.0 ** 2)


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
				print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
		else:
			total_downloaded = blocknum * block_size
			status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

			if total_size > 0:
				percent_downloaded = total_downloaded * 100.0 / total_size
				# use carriage return plus sys.stdout to overwrite stdout
				try:
					download_rate = total_downloaded / (time.time() - self.start_time)
					estimated_time = (total_size - total_downloaded) / download_rate
				except ZeroDivisionError:
					download_rate = 0
					estimated_time = 0
				estimated_minutes = int(estimated_time / 60.0)
				estimated_seconds = estimated_time - estimated_minutes * 60.0
				status += "{:3.2f}".format(percent_downloaded) + " %  " + \
				          "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
				          "{:2.0f}".format(estimated_minutes) + " min " + \
				          "{:2.0f}".format(estimated_seconds) + " sec "
			status += "        \r"
			sys.stdout.write(status)


def download(url, download_file):
	"""
	Download a file from a url
	"""

	try:
		print("Downloading " + url)
		file, headers = urlretrieve(url, download_file, reporthook=ReportHook().report)
		# print final return to start new line of stdout
		print("\n")
	except EnvironmentError:
		print("WARNING: Unable to download " + url)


def download_unpack_tar(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		tarfile_handle = tarfile.open(download_file)
		tarfile_handle.extractall(path=folder)
		tarfile_handle.close()
	except (EnvironmentError, tarfile.ReadError):
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)


def download_unpack_zip(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		zipfile_handle = zipfile.ZipFile(download_file)
		zipfile_handle.extractall(path=folder)
		zipfile_handle.close()
	except EnvironmentError:
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)


setuptools.setup(
	name="metawibele",
	author=AUTHOR,
	author_email=MAINTAINER_EMAIL,
	version=VERSION,
	license="MIT",
	description="MetaWIBELE: Workflow to Identify novel Bioactive Elements in microbiome",
	long_description="MetaWIBELE (Workflow to Identify novel Bioactive Elements in the microbiome) " +
	                 "is a computational pipeline that identifies novel bioactive microbial gene products " +
					 "from metagenomes and finds new immunomodulatory gene families, " +
					 "especially targeting secreted/extracellular proteins to enrich for likely host interactors. " +
					 "The prioritized list of gene products can be further used for downstream experimental validation.",
	url="https://github.com/biobakery/metawibele",
	keywords=['microbial', 'microbiome', 'bioinformatics', 'microbiology', 'metagenomic', 'metatranscriptomic',
	         'gene catalog', 'bioactive', 'assembly',
	          'metawibele', 'anadama2'],
	platforms=['Linux', 'MacOS'],
	classifiers=[
		"Programming Language :: Python",
		"Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Operating System :: MacOS",
		"Operating System :: Unix",
		"Programming Language :: Python :: 3.6",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	#install_requires=['anadama2>=0.7.4'],
	packages=setuptools.find_packages(),
	#cmdclass={'install': Install},
	entry_points={
		'console_scripts': [
			'metawibele = metawibele.metawibele:main',
			'metawibele_characterize = metawibele.workflows.characterize:main',
			'metawibele_prioritize = metawibele.workflows.prioritize:main',
			'metawibele_preprocess = metawibele.workflows.preprocess:main',
			'metawibele_abundance_annotator = metawibele.characterize.abundance_annotator:main',
			'metawibele_antiSMASH_annotator = metawibele.characterize.antiSMASH_annotator:main',
			'metawibele_ddi_DOMINE_ExpAtlas = metawibele.characterize.ddi_DOMINE_ExpAtlas:main',
			'metawibele_ddi_DOMINE_SIFTS = metawibele.characterize.ddi_DOMINE_SIFTS:main',
			'metawibele_ddi_DOMINE_ann = metawibele.characterize.ddi_DOMINE_ann:main',
			'metawibele_ddi_DOMINE_protein = metawibele.characterize.ddi_DOMINE_protein:main',
			'metawibele_ddi_DOMINE_protein_family = metawibele.characterize.ddi_DOMINE_protein_family:main',
			'metawibele_denovo_TM_SP = metawibele.characterize.denovo_TM_SP:main',
			'metawibele_finalize_annotation = metawibele.characterize.finalize_annotation:main',
			'metawibele_interproscan_annotator = metawibele.characterize.interproscan_annotator:main',
			'metawibele_interproscan_pfam_protein_family = metawibele.characterize.interproscan_pfam_protein_family:main',
			'metawibele_interproscan_phobius_protein_family = metawibele.characterize.interproscan_phobius_protein_family:main',
			'metawibele_interproscan_protein = metawibele.characterize.interproscan_protein:main',
			'metawibele_interproscan_protein_family = metawibele.characterize.interproscan_protein_family:main',
			'metawibele_interproscan_signalp_protein_family = metawibele.characterize.interproscan_signalp_protein_family:main',
			'metawibele_interproscan_tmhmm_protein_family = metawibele.characterize.interproscan_tmhmm_protein_family:main',
			'metawibele_maaslin2 = metawibele.characterize.maaslin2:main',
			'metawibele_maaslin2_annotator = metawibele.characterize.maaslin2_annotator:main',
			'metawibele_maaslin2_collection = metawibele.characterize.maaslin2_collection:main',
			'metawibele_maaslin2_summary = metawibele.characterize.maaslin2_summary:main',
			'metawibele_metadata_format = metawibele.characterize.metadata_format:main',
			'metawibele_msp_protein_family = metawibele.characterize.msp_protein_family:main',
			'metawibele_mspminer_msp = metawibele.characterize.mspminer_msp:main',
			'metawibele_mspminer_msp_taxonomy_annotation = metawibele.characterize.mspminer_msp_taxonomy_annotation:main',
			'metawibele_mspminer_msp_uniref_annotation = metawibele.characterize.mspminer_msp_uniref_annotation:main',
			'metawibele_mspminer_protein = metawibele.characterize.mspminer_protein:main',
			'metawibele_mspminer_protein_family = metawibele.characterize.mspminer_protein_family:main',
			'metawibele_mspminer_protein_family_taxonomy = metawibele.characterize.mspminer_protein_family_taxonomy:main',
			'metawibele_pfam2go = metawibele.characterize.pfam2go:main',
			'metawibele_psortb_annotator = metawibele.characterize.psortb_annotator:main',
			'metawibele_psortb_protein = metawibele.characterize.psortb_protein:main',
			'metawibele_psortb_protein_family = metawibele.characterize.psortb_protein_family:main',
			'metawibele_sum_to_protein_family_abundance = metawibele.characterize.sum_to_protein_family_abundance:main',
			'metawibele_sum_to_protein_family_stratified_abundance = metawibele.characterize.sum_to_protein_family_stratified_abundance:main',
			'metawibele_summary_all_annotation = metawibele.characterize.summary_all_annotation:main',
			'metawibele_summary_function_annotation = metawibele.characterize.summary_function_annotation:main',
			'metawibele_summary_protein_family_uniref_annotation = metawibele.characterize.summary_protein_family_uniref_annotation:main',
			'metawibele_summary_protein_uniref_annotation = metawibele.characterize.summary_protein_uniref_annotation:main',
			'metawibele_uniref_annotator = metawibele.characterize.uniref_annotator:main',
			'metawibele_uniref_annotator_stat = metawibele.characterize.uniref_annotator_stat:main',
			'metawibele_uniref_protein = metawibele.characterize.uniref_protein:main',
			'metawibele_uniref_protein_family = metawibele.characterize.uniref_protein_family:main',
			'metawibele_cluster_prioritization = metawibele.prioritize.cluster_prioritization:main',
			'metawibele_filter_prioritization = metawibele.prioritize.filter_prioritization:main',
			'metawibele_finalize_prioritization = metawibele.prioritize.finalize_prioritization:main',
			'metawibele_quantify_prioritization = metawibele.prioritize.quantify_prioritization:main',
			'metawibele_combine_gene_sequences = metawibele.tools.combine_gene_sequences:main',
			'metawibele_extract_cluster = metawibele.tools.extract_cluster_CD_hit:main',
			'metawibele_extract_complete_ORF_seq = metawibele.tools.extract_complete_ORF_seq:main',
			'metawibele_extract_non_redundance_seq = metawibele.tools.extract_non_redundance_seq:main',
			'metawibele_extract_protein_coding_genes = metawibele.tools.extract_protein_coding_genes:main',
			'metawibele_format_contig_sequences = metawibele.tools.format_contig_sequences:main',
			'metawibele_format_protein_sequences = metawibele.tools.format_protein_sequences:main',
			'metawibele_gene_abundance = metawibele.tools.gene_abundance:main',
			'metawibele_gene_abundance_indexRef = metawibele.tools.gene_abundance_indexRef:main',
			'metawibele_gene_catalog_abundance = metawibele.tools.gene_catalog_abundance:main',
			'metawibele_abundance_RPK = metawibele.common.abundance_RPK:main',
			'metawibele_abundance_RPK_gene = metawibele.common.abundance_RPK_gene:main',
			'metawibele_abundance_filtering = metawibele.common.abundance_filtering:main',
			'metawibele_abundance_normalization = metawibele.common.abundance_normalization:main',
			'metawibele_abundance_smoothing = metawibele.common.abundance_smoothing:main',
			'metawibele_combine_abundance_annotation = metawibele.common.combine_abundance_annotation:main',
			'metawibele_extract_abundance_feature_subset = metawibele.common.extract_abundance_feature_subset:main',
			'metawibele_extract_abundance_sample_subset = metawibele.common.extract_abundance_sample_subset:main',
			'metawibele_filter_clusters = metawibele.common.filter_clusters:main',
			'metawibele_filter_prevalence = metawibele.common.filter_prevalence:main',
			'metawibele_join_family_abundance = metawibele.common.join_family_abundance:main',
			'metawibele_split_family_abundance = metawibele.common.split_family_abundance:main',
			'metawibele_transpose = metawibele.common.transpose:main',
			'metawibele_split_fasta_file = metawibele.common.split_fasta_file:main',
			'metawibele_download_config = metawibele.common.download_config_file:main',
			'metawibele_download_database = metawibele.common.download_database:main',
			'metawibele_prepare_uniprot_taxonomy = metawibele.common.prepare_uniprot_taxonomy:main',
			'metawibele_prepare_uniprot_annotation = metawibele.common.prepare_uniprot_annotation:main',
			'metawibele_prepare_uniref_annotation = metawibele.common.prepare_uniref_annotation:main',
			'metawibele_prepare_domain_databases = metawibele.common.prepare_domain_databases:main',
			'metawibele_extract_uniref_maps = metawibele.common.extract_uniref_maps:main',
			'metawibele_check_install = metawibele.common.check_install:main'
		]},
	package_data={
		'metawibele': [
			'workflows/*.py',
			'metawibele.cfg',
			'configs/*',
			'data/domain/*',
			'data/uniref/*',
			'data/misc/*',
			'Rscripts/*'
		]},
	#data_files = [
	#		("metawibele/examples/", glob("examples/*"))
	#	],
	scripts=glob('metawibele/workflows/*py'),
	zip_safe=False
)

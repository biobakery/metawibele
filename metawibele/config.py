#!/usr/bin/env python

"""
MetaWIBELE: config module
Configuration settings

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

# try to import the python2 ConfigParser
# if unable to import, then try to import the python3 configparser
try:
    import ConfigParser as configparser
except ImportError:
    import configparser

import logging


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Config metawibele
"""

def get_args (): 
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-c', "--config",
                        help='input global config file',
                        required=False,
						default="none")
    values = parser.parse_args()
    return values
# get_args


##default option for MetaWIBELE
version = '0.0.1'
log_level = 'DEBUG'
verbose = 'DEBUG'
PROTEIN_FAMILY_ID = "familyID"
PROTEIN_ID = "seqID"
#arg_values = get_args ()

# name global logging instance
logger=logging.getLogger(__name__)

def log_settings():
	"""
	Write to the log file the config settings for the run
	"""

	lines = []
	lines.append("DATABASE SETTINGS")
	lines.append("nucleotide database folder = " + nucleotide_database)
	lines.append("protein database folder = " + protein_database)
	if pathways_database_part1:
		lines.append("pathways database file 1 = " + pathways_database_part1)
		lines.append("pathways database file 2 = " + pathways_database_part2)
	else:
		lines.append("pathways database file = " + pathways_database_part2)
	lines.append("utility mapping database folder = " + utility_mapping_database)
	lines.append("")

	lines.append("RUN MODES")
	lines.append("resume = " + str(resume))
	lines.append("verbose = " + str(verbose))
	lines.append("bypass prescreen = " + str(bypass_prescreen))
	lines.append("bypass nucleotide index = " + str(bypass_nucleotide_index))
	lines.append("bypass nucleotide search = " + str(bypass_nucleotide_search))
	lines.append("bypass translated search = " + str(bypass_translated_search))
	lines.append("translated search = " + translated_alignment_selected)
	lines.append("pick frames = " + pick_frames_toggle)
	lines.append("threads = " + str(threads))
	lines.append("")

	lines.append("SEARCH MODE")
	lines.append("search mode = " + search_mode)
	lines.append("identity threshold = " + str(identity_threshold))
	lines.append("")

	lines.append("ALIGNMENT SETTINGS")
	lines.append("evalue threshold = " + str(evalue_threshold))
	lines.append("prescreen threshold = " + str(prescreen_threshold))
	lines.append("translated subject coverage threshold = " + str(translated_subject_coverage_threshold))
	lines.append("translated query coverage threshold = " + str(translated_query_coverage_threshold))
	lines.append("")

	lines.append("PATHWAYS SETTINGS")
	lines.append("minpath = " + minpath_toggle)
	lines.append("xipe = " + xipe_toggle)
	lines.append("gap fill = " + gap_fill_toggle)
	lines.append("")

	lines.append("INPUT AND OUTPUT FORMATS")
	lines.append("input file format = " + input_format)
	lines.append("output file format = " + output_format)
	lines.append("output max decimals = " + str(output_max_decimals))
	lines.append("remove stratified output = " + str(remove_stratified_output))
	lines.append("remove column description output = " + str(remove_column_description_output))
	lines.append("log level = " + log_level)
	lines.append("")

	logger.info("\nRun config settings: \n\n" + "\n".join(lines))


def update_user_edit_config_file_single_item(section, name, value):
	"""
	Update the settings to the user editable config file for one item
	"""

	new_config_items = {section: {name: value}}

	update_user_edit_config_file(new_config_items)

	print("MetaWIBELE configuration file updated: " + section + " : " + name + " = " + str(value))


def update_user_edit_config_file(new_config_items):
	"""
	Update the settings to the user editable config file
	"""

	config = configparser.RawConfigParser()

	# start with the current config settings
	config_items = read_user_edit_config_file(full_path_user_edit_config_file)

	# update with the new config items
	for section in new_config_items:
		for name, value in new_config_items[section].items():
			if section in config_items:
				if name in config_items[section]:
					config_items[section][name] = value
				else:
					sys.exit("ERROR: Unable to add new name ( " + name +
					         " ) to existing section ( " + section + " ) to " +
					         " config file: " + full_path_user_edit_config_file)
			else:
				sys.exit("ERROR: Unable to add new section ( " + section +
				         " ) to config file: " + full_path_user_edit_config_file)

	for section in config_items:
		config.add_section(section)
		for name, value in config_items[section].items():
			value = str(value)
			if "file" in section or "folder" in section:
				# convert to absolute path if needed
				if not os.path.isabs(value):
					value = os.path.abspath(value)
			config.set(section, name, value)

	try:
		file_handle = open(full_path_user_edit_config_file, "wt")
		config.write(file_handle)
		file_handle.close()
	except EnvironmentError:
		sys.exit("Unable to write to the MetaWIBELE config file.")


def read_user_edit_config_file(full_path_user_edit_config_file):
	"""
	Read the settings from the config file
	"""

	config = configparser.ConfigParser()

	try:
		config.read(full_path_user_edit_config_file)
	except EnvironmentError:
		sys.exit("Unable to read from the config file: " + full_path_user_edit_config_file)

	# read through all of the sections
	config_items = {}
	for section in config.sections():
		config_list = config.items(section)
		config_items[section] = {}
		for name, value in config_list:
			if "file" in section or "folder" in section:
				# if not absolute path, then return absolute path relative to this folder
				if not os.path.isabs(value):
					value = os.path.abspath(os.path.join(os.path.dirname(full_path_user_edit_config_file), value))
			config_items[section][name] = value

	return config_items


def get_item(config_items, section, name, type=None):
	"""
	Get the item from the dictionary of section/names from the user edit config file
	"""

	# try to obtain the value from the config dictionary
	try:
		value = config_items[section][name]
	except KeyError:
		sys.exit("CRITICAL ERROR: Unable to load value from " + full_path_user_edit_config_file +
		         " . \nItem not found. \nItem should be in section (" + section + ") with name (" + name + ").")

	# if present, try to change the value type
	if type:
		try:
			if type == "string":
				value = str(value)
			elif type == "int":
				value = int(value)
			elif type == "float":
				value = float(value)
			elif type == "bool":
				if value in ["False", "false", "F", "f"]:
					value = False
				elif value in ["True", "true", "T", "t"]:
					value = True
				else:
					raise ValueError
		except ValueError:
			sys.exit("CRITICAL ERROR: Unable to load value from " + full_path_user_edit_config_file +
			         " . \nItem found in section (" + section + ") with name (" + name + "). " +
			         "\nItem is not of type (" + type + ").")

	return value


# User config file
#metawibele_install_directory = "/n/home00/yancong/projects/R24_HMBR/src/assembly-based/metawibele/" 
metawibele_install_directory = os.path.dirname(os.path.abspath(__file__))

user_edit_config_file = os.getcwd() + "/metawibele.cfg"
if not os.path.exists (user_edit_config_file):
	user_edit_config_file = metawibele_install_directory + "/metawibele.cfg"
#full_path_user_edit_config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
#												user_edit_config_file)
full_path_user_edit_config_file = user_edit_config_file

# get the base settings from the user edit config file
config_items = read_user_edit_config_file(full_path_user_edit_config_file)

## databases ##
database_directory = metawibele_install_directory + "/data/"
taxonomy_database_choices = ["uniprot_taxonomy.map.tsv","uniprot_taxaID_bac-arc-vir.tsv","uniprot_taxaID_mammalia.tsv"]
taxonomy_database = database_directory + taxonomy_database_choices[0]
microbiome_taxa = database_directory + taxonomy_database_choices[1]
mammalia_taxa = database_directory + taxonomy_database_choices[2]
pdb_database_choices = ["pdb_chain_taxonomy.tsv","pdb_chain_pfam.tsv"]
pdb_taxonomy = database_directory + pdb_database_choices[0]
pdb_pfam = database_directory + pdb_database_choices[1]
uniref_database_choices = ['uniref90.ann.all.tsv', 'uniref50.ann.all.tsv', 'map_UniRef90_UniRef50.dat', 'uniref90.fasta.dmnd', 'uniref50.fasta.dmnd', 'uniref90.fasta', 'uniref50.fasta']
uniref_database = database_directory + uniref_database_choices[0]
uniref50_database = database_directory + uniref_database_choices[1]
uniref_map = database_directory + uniref_database_choices[2]
uniref_dmnd = database_directory + uniref_database_choices[3]
uniref50_dmnd = database_directory + uniref_database_choices[4]
pfam_database_choices = ["Pfam_ann.tsv","uniprot_human_pfam.tsv"]
pfam_database = database_directory + pfam_database_choices[0]
human_pfam_database = database_directory + pfam_database_choices[1]
pfam2go_database_choices = ['Pfam2GO.txt']
pfam2go_database = database_directory + pfam2go_database_choices[0]
go_database_choices = ['go.obo']
go_database = database_directory + go_database_choices[0]
interaction_database_choices = ['INTERACTION.txt']
interaction_database = database_directory + interaction_database_choices[0]
Expression_Atlas_database = [database_directory + "32_Uhlen_Lab_colon_rectum.tsv", 
							database_directory + "Encode_sigmoid_colon.tsv",
							database_directory + "FANTOM5_colon_rectum.tsv",
							database_directory + "GTEx_sigmoid_transverse_colon.tsv",
							database_directory + "Human_protein_Atlas_colon_rectum.tsv",
							database_directory + "Human_proteome_map_colon_rectum.tsv",
							database_directory + "Illumina_Body_Map_colon.tsv"]
vignettes_database = database_directory + "vignettes_proteins.tsv"
antiSMASH_database = database_directory + "BGC_genes_unirefID.tsv"
#essentialGene_database = database_directory + "DEG_genes_unirefID.tsv"
essentialGene_database = database_directory + "essential_genes_unirefID.tsv"

## config files and computing resources
config_directory = metawibele_install_directory + "/config/"
threads = get_item (config_items, "computation", "threads", "int")

## Input and output ##
# input folder and files
basename = get_item(config_items, "output", "basename", "string")
working_dir = get_item(config_items, "output", "working_dir", "string")
gene_catalog_dir = working_dir + "/input/"
output_dir = working_dir + "/output/" 
annotation_dir = working_dir + "/" + "characterization"
cluster_dir = annotation_dir + "/" + "cluster"
homology_dir = annotation_dir + "/" + "homology_annotation"
sequence_dir = annotation_dir + "/" + "sequence_annotation"
abundance_dir = annotation_dir + "/" + "abundance_annotation"
priority_dir = working_dir + "/" + "prioritization"

study = get_item(config_items, "input_dataset", "study", "string")
metadata = get_item(config_items, "input_dataset", "metadata", "string")
rna_metadata = get_item(config_items, "input_dataset", "rna_metadata", "string")
sample_list = get_item(config_items, "input_dataset", "sample_list", "string")
protein_seq = get_item(config_items, "input_dataset", "protein_seq", "string")
gene_catalog = get_item(config_items, "input_dataset", "gene_catalog", "string")
gene_catalog_nuc = get_item(config_items, "input_dataset", "gene_catalog_nuc", "string")
gene_catalog_prot = get_item(config_items, "input_dataset", "gene_catalog_prot", "string")
gene_catalog_count = get_item(config_items, "input_dataset", "gene_catalog_count", "string")

protein_family = cluster_dir + "/" + basename + "_proteinfamilies.clstr"
protein_family_seq = cluster_dir + "/" + basename + "_proteinfamilies.centroid.fasta"
protein_family_relab = abundance_dir + "/" + basename + "_proteinfamilies_relab.tsv"
protein_family_ann = annotation_dir + "/" + basename + "_proteinfamilies_annotation.tsv" 
protein_family_attr = annotation_dir + "/" + basename + "_proteinfamilies_annotation.attribute.tsv"
unsupervised_rank = priority_dir + "/" + basename + "_unsupervised_prioritization.rank.table.tsv"
supervised_rank = priority_dir + "/" + basename + "_supervised_prioritization.rank.table.tsv"
unsupervised_priority = priority_dir + "/" + basename + "_unsupervised_prioritization.priority.table.tsv"
supervised_priority = priority_dir + "/" + basename + "_supervised_prioritization.priority.table.tsv"


## characterization
#tshld_consistency = 0.75	# the minimum annotation consistency in one protein family
tshld_consistency = get_item(config_items, "characterization", "tshld_consistency", "float") # the minimum annotation consistency in one protein family

# CD-hit
cd_hit_memory = get_item(config_items, "CD-hit", "cd_hit_memory", "int")
cd_hit_protein_opts = get_item(config_items, "CD-hit", "cd_hit_protein_opts", "string")  # used for protein families
cd_hit_gene_opts = get_item(config_items, "CD-hit", "cd_hit_gene_opts", "string") # ued for gene catalogs
featureCounts_opts = get_item(config_items, "CD-hit", "featurecounts_opts", "string") + " -T " + str(threads)

# diamond options
diamond_database_extension = get_item(config_items, "diamond", "diamond_database_extension", "string")
diamond_cmmd_protein_search = get_item(config_items, "diamond", "diamond_cmmd_protein_search", "string")
diamond_cmmd_nucleotide_search = get_item(config_items, "diamond", "diamond_cmmd_nucleotide_search", "string")
diamond_version={
    "flag" : "--version",
    "major" : 0,
    "minor" : 8,
    "second minor" : 22,
    "line" : 0,
    "column" : 2}
diamond_identity = get_item(config_items, "diamond", "diamond_identity", "float")
diamond_query_coverage = get_item(config_items, "diamond", "diamond_query_coverage", "float")
diamond_mutual_coverage = get_item(config_items, "diamond", "diamond_mutual_coverage", "float")

# homology-based
tshld_identity = get_item(config_items, "homology-based", "tshld_identity", "float")    # the minimum identity of homology
tshld_coverage = get_item(config_items, "homology-based", "tshld_coverage", "float") 	# the minimum coverage of homology
taxa_source = get_item(config_items, "homology-based", "taxa_source", "string") 		# the source of taxa for one protein family, representatives vs. LCA

# interporscan
interproscan_appl = get_item(config_items, "interproscan", "interproscan_appl", "string")
interproscan_type = []
tmp = interproscan_appl.split(",")
for item in tmp:
	if item == "Phobius":
		item = "phobius.signaling"
		item = "interpro." + item + ".tsv"
		interproscan_type.append(item)
		item = "phobius.transmembrane"
		item = "interpro." + item + ".tsv"
		interproscan_type.append(item)
		continue
	if item == "Pfam":
		item = "PfamDomain"
	if item == "SignalP":
		item = "signalp.signaling"
	if item == "TMHMM":
		item = "tmhmm.transmembrane"
	item = "interpro." + item + ".tsv"
	interproscan_type.append(item)


# MSP
#msp_dir = get_item(config_items, "output_folders", "msp_dir", "string")
#msp_dir = abundance_dir + "/MSPminer/"
#msp_config =  get_item(config_items, "MSP", "msp_config", "string")
tshld_unclassified = get_item(config_items, "MSP", "tshld_unclassified", "float")	# the minimum percentile of unclassified in one MSP
tshld_diff = get_item(config_items, "MSP", "tshld_diff", "float")			        # the minimum difference between most and second dominant taxa in the MSP
tshld_lca = get_item(config_items, "MSP", "tshld_lca", "float")				        # the minimum percentile of consistent taxon for protein family to get LCA
taxa_final = get_item(config_items, "MSP", "taxa_final", "string")				    # the source of taxa for one protein family, representatives vs. LCA

# abundance
normalization = get_item(config_items, "abundance", "normalization", "string") 		# the method for abundance normalization 

# maaslin2 options
maaslin2_dir = get_item(config_items, "maaslin2", "maaslin2_output", "string")
maaslin2_dir = abundance_dir + "/" + maaslin2_dir + "/maaslin2_output/" 
#maaslin2_dir = "/n/scratchlfs/huttenhower_lab/yancong/assembly-based/HMP2/combined/HUMAnN2/maaslin2_output/"
contrast_status = get_item(config_items, "maaslin2", "contrast_status", "string")
tmp = contrast_status.split("|")
contrast_status = {}
for item in tmp:
	tmp1 = item.split(":")
	contrast_status[tmp1[0]] = tmp1[1]
ref_status = get_item(config_items, "maaslin2", "ref_status", "string")
tmp = ref_status.split("|")
ref_status = {}
for item in tmp:
	tmp1 = item.split(":")
	ref_status[tmp1[0]] = tmp1[1]
#fixed_effect = "diagnosis,age,antibiotic,immunosuppressant,mesalamine,steroids"
fixed_effects = get_item(config_items, "maaslin2", "fixed_effects", "string")
random_effects = get_item(config_items, "maaslin2", "random_effects", "string")
nested_effects = get_item(config_items, "maaslin2", "nested_effects", "string")
#abundance_detection_level = 0	# the detectable level of abundance
abundance_detection_level = get_item(config_items, "maaslin2", "abundance_detection_level", "float")
#tshld_prevalence = 0.10		# the minimum prevalence of protein family abundance across samples
tshld_prevalence = get_item(config_items, "maaslin2", "tshld_prevalence", "float")
#tshld_qvalue = 0.05			# the maximum q-value for significant results
tshld_qvalue = get_item(config_items, "maaslin2", "tshld_qvalue", "float")
#effect_size = 'coef'			# the source of effect size, coef vs. Cohen's d
effect_size = get_item(config_items, "maaslin2", "effect_size", "string")
#maaslin2_cores = 25			# the number of threads used by maaslin2
maaslin2_cores = get_item(config_items, "maaslin2", "maaslin2_cores", "int")
#maaslin2_cmmd = "~/usr/Maaslin2/R/Maaslin2.R"
maaslin2_cmmd = get_item(config_items, "maaslin2", "maaslin2_cmmd", "string")
#maaslin2_utils = get_item(config_items, "maaslin2", "maaslin2_utils", "string")
maaslin2_utils = metawibele_install_directory + "/characterize/maaslin2_utils.r"
#pcl_utils = get_item(config_items, "maaslin2", "pcl_utils", "string")
pcl_utils = metawibele_install_directory + "/characterize/pcl_utils.r"
#transpose_cmmd = get_item(config_items, "maaslin2", "transpose_cmmd", "string")
transpose_cmmd =  metawibele_install_directory + "/common/transpose.py"
min_abundance = get_item(config_items, "maaslin2", "min_abundance", "float")
min_prevalence = get_item(config_items, "maaslin2", "min_prevalence", "float")
max_significance = get_item(config_items, "maaslin2", "max_significance", "float")
normalization = get_item(config_items, "maaslin2", "normalization", "string")
transform = get_item(config_items, "maaslin2", "transform", "string")
analysis_method = get_item(config_items, "maaslin2", "analysis_method", "string")
plot_heatmap = get_item(config_items, "maaslin2", "plot_heatmap", "string")
plot_scatter = get_item(config_items, "maaslin2", "plot_scatter", "string")
maaslin2_cmmd_opts=["--min_abundance", min_abundance, "--min_prevalence", min_prevalence, "--max_significance", max_significance, "--normalization", normalization,  "--transform", transform, "--analysis_method", analysis_method, "--cores", maaslin2_cores, "--fixed_effects", fixed_effects, "--random_effects", random_effects, "--plot_heatmap", plot_heatmap, "--plot_scatter", plot_scatter]


## prioritization thresholds 
#tshld_abund = 75
#tshld_abund = 75
#tshld_prev = 75
#tshld_priority = 0.25		# the top percentile for priority
tshld_priority = get_item(config_items, "prioritization", "tshld_priority", "float")
#tshld_score = 0.5			# the minimum score for priority 
tshld_priority_score = get_item(config_items, "prioritization", "tshld_priority_score", "float")


## assembly workflow

# bowtie2 options and threshold
bowtie2_large_index_threshold=4000000000
bowtie2_index_ext_list=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",
    ".rev.1.bt2",".rev.2.bt2"]
bowtie2_large_index_ext=".1.bt2l"
bowtie2_version={
    "flag" : "--version",
    "major" : 2,
    "minor" : 2,
    "line" : 0,
    "column" : 2}
bowtie2_align_opts=["--threads", threads] # "--threads "+str(threads)
bowtie2_index_name="_bowtie2_index"


if __name__=='__main__':
	pass
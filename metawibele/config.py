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
import re

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


## default option for MetaWIBELE ##
version = '0.3.2'
log_level = 'DEBUG'
verbose = 'DEBUG'

# name global logging instance
logger = logging.getLogger(__name__)


## constant values ##
PROTEIN_FAMILY_ID = "familyID"
PROTEIN_ID = "seqID"
c_metedata_delim = "."	 # nested metadata, e.g. CD.dysbiosis
c_strat_delim = "|" 	 # strantified item, e.g. Cluster_1000010|Bacteroides dorei
c_taxon_delim = "|"	 	 # taxonomic lineage, e.g. g__Faecalibacterium|s__Faecalibacterium_prausnitzii|t__Faecalibacterium_prausnitzii_A2-165
c_multiname_delim = ";"	 # multiple ietms, e.g. PF00482;PF01841 
c_msp_unknown = "msp_unknown"


## User config file ##
metawibele_install_directory = os.path.dirname(os.path.abspath(__file__))
config_directory = os.path.join(metawibele_install_directory, "configs")

user_edit_config_file = os.path.join(os.getcwd(), "metawibele.cfg")
if not os.path.exists (user_edit_config_file):
	#user_edit_config_file = metawibele_install_directory + "/metawibele.cfg"
	user_edit_config_file = os.path.join(config_directory, "metawibele.cfg")
full_path_user_edit_config_file = user_edit_config_file

# get the base settings from the user edit config file
config_items = read_user_edit_config_file(full_path_user_edit_config_file)


## Databases ##
# installed data
database_directory = os.path.join(metawibele_install_directory, "data")
taxonomy_directory = os.path.join(database_directory, "taxonomy")
taxonomy_database_choices = ["uniprot_taxonomy.txt.gz","uniprot_taxaID_bac-arc-vir.txt.gz","uniprot_taxaID_mammalia.txt.gz"]
taxonomy_database = os.path.join(taxonomy_directory,  taxonomy_database_choices[0])
microbiome_taxa = os.path.join(taxonomy_directory, taxonomy_database_choices[1])
mammalia_taxa = os.path.join(taxonomy_directory, taxonomy_database_choices[2])

domain_directory = os.path.join(database_directory, "domain")
pdb_database_choices = ["pdb_chain_taxonomy.txt.gz","pdb_chain_pfam.txt.gz"]
pdb_taxonomy = os.path.join(domain_directory, pdb_database_choices[0])
pdb_pfam = os.path.join(domain_directory, pdb_database_choices[1])
pfam_database_choices = ["pfam_descriptions.txt.gz","uniprot_human_pfam.txt.gz"]
pfam_database = os.path.join(domain_directory, pfam_database_choices[0])
human_pfam_database = os.path.join(domain_directory, pfam_database_choices[1])
pfam2go_database_choices = ['pfam2GO.txt.gz']
pfam2go_database = os.path.join(domain_directory, pfam2go_database_choices[0])
interaction_database_choices = ['INTERACTION.txt.gz']
interaction_database = os.path.join(domain_directory, interaction_database_choices[0])
Expression_Atlas_database = [os.path.join(domain_directory, "32_Uhlen_Lab_colon_rectum.txt.gz"), 
							os.path.join(domain_directory, "Encode_sigmoid_colon.txt.gz"),
							os.path.join(domain_directory, "FANTOM5_colon_rectum.txt.gz"),
							os.path.join(domain_directory, "GTEx_sigmoid_transverse_colon.txt.gz"),
							os.path.join(domain_directory,"Human_protein_Atlas_colon_rectum.txt.gz"),
							os.path.join(domain_directory, "Human_proteome_map_colon_rectum.txt.gz"),
							os.path.join(domain_directory, "Illumina_Body_Map_colon.txt.gz")]

misc_directory = os.path.join(database_directory, "misc")
vignettes_database = os.path.join(misc_directory, "vignette_function.tsv")

# uniref databases
uniref_database_dir = get_item (config_items, "uniref", "database", "string")
if uniref_database_dir.lower() == "none" or uniref_database_dir == "":
	uniref_database_dir = os.path.join(database_directory, "uniref")
uniref_database_choices = ['uniref90.fasta.dmnd', 'uniref90.ann.tsv.gz']
uniref_dmnd = os.path.join(uniref_database_dir, uniref_database_choices[0])
'''
files_path1 = []
if os.isdir(os.path.join(database_directory, "uniref")):
	path = os.path.join(database_directory, "uniref")
	files_path1 = [os.path.abspath(x) for x in os.listdir(path)]
files_path2 = []
if os.isdir(get_item (config_items, "uniref", "database", "string")):
	path = get_item (config_items, "uniref", "database", "string")
	files_path2 = [os.path.abspath(x) for x in os.listdir(path)]
uniref_database = files_path1 + files_path2
'''
uniref_database = [os.path.join(uniref_database_dir, "map_Protein_names_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_Gene_names_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_UniProtKB_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_Tax_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_TaxID_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_Rep_Tax_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_Rep_TaxID_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_go_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_ko_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_eggnog_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_level4ec_uniref90.txt.gz"),
                        os.path.join(uniref_database_dir, "map_pfam_uniref90.txt.gz")]


## Computing resources ##
threads = get_item (config_items, "computation", "threads", "int")
memory = get_item(config_items, "computation", "memory", "int")
time = get_item(config_items, "computation", "time", "int")


## Input and output ##
# input folder and files
basename = get_item(config_items, "output", "basename", "string")
working_dir = get_item(config_items, "output", "output_dir", "string")
if working_dir.lower() == "none" or working_dir.lower() == "":
	working_dir = os.getcwd()
annotation_dir = os.path.join(working_dir, "characterization")
cluster_dir = os.path.join(annotation_dir, "clustering")
global_homology_dir = os.path.join(annotation_dir, "global_homology_annotation")
domain_motif_dir = os.path.join(annotation_dir, "doamin_motif_annotation")
abundance_dir = os.path.join(annotation_dir, "abundance_annotation")
priority_dir = os.path.join(working_dir, "prioritization")

study = get_item(config_items, "input", "study", "string")
metadata = get_item(config_items, "input", "metadata", "string")
sample_list = get_item(config_items, "input", "sample_list", "string")
gene_catalog = get_item(config_items, "input", "gene_catalog", "string")
gene_catalog_prot = get_item(config_items, "input", "gene_catalog_prot", "string")
gene_catalog_count = get_item(config_items, "input", "gene_catalog_count", "string")

protein_family = os.path.join(annotation_dir, basename + "_proteinfamilies.clstr")
protein_family_prot_seq = os.path.join(annotation_dir, basename + "_proteinfamilies.centroid.faa")
protein_family_nuc_seq = os.path.join(annotation_dir, basename + "_proteinfamilies.centroid.fna")
protein_family_relab = os.path.join(annotation_dir, basename + "_proteinfamilies_relab.tsv")
protein_family_ann = os.path.join(annotation_dir, basename + "_proteinfamilies_annotation.tsv") 
protein_family_attr = os.path.join(annotation_dir, basename + "_proteinfamilies_annotation.attribute.tsv")
unsupervised_rank = os.path.join(priority_dir, basename + "_unsupervised_prioritization.rank.table.tsv")
supervised_rank =  os.path.join(priority_dir, basename + "_supervised_prioritization.rank.table.tsv")


## characterization ##
tshld_consistency = 0.75	# the minimum annotation consistency in one protein family

# CD-hit
# memory used for CD-hit
cd_hit_memory = memory 
cd_hit_prot_opts = "-d 100 -c 0.9 -aL 0.8 -aS 0.8 -G 0 -M 0 -B 0"	# clustering protein families
cd_hit_gene_opts = "-d 100 -c 0.95 -aS 0.8 -G 0 -M 0 -B 0"
featureCounts_opts = " -F SAF "

# diamond options
diamond_database_extension = ".dmnd"
diamond_cmmd_protein_search = "blastp"
diamond_cmmd_nucleotide_search = "blastx"
diamond_identity = 0.9			# identity threshold for uniref90 strong homologies
diamond_query_coverage = 0.8 	# query and mutual coverage threshold of uniref90 strong homologies
diamond_mutual_coverage = 0.8
diamond_version={
    "flag" : "--version",
    "major" : 0,
    "minor" : 8,
    "second minor" : 22,
    "line" : 0,
    "column" : 2}

# protein family
tshld_identity = 0.25	# the minimum identity of homology
tshld_coverage = 0.25	# the minimum coverage of homology
taxa_source = "Rep"		# the source of taxa for one protein family, representatives vs. LCA

# interporscan
interproscan_cmmd = get_item(config_items, "interproscan", "interproscan_cmmd", "string")
interproscan_cmmd = re.sub("\"", "", interproscan_cmmd)
interproscan_appl = get_item(config_items, "interproscan", "interproscan_appl", "string")
if interproscan_appl.lower() == "none" or interproscan_appl.lower() == "":
	interproscan_appl = "CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM,Phobius,SignalP,TMHMM"
interproscan_appl = re.sub("\"", "", interproscan_appl)
split_number = get_item(config_items, "interproscan", "split_number", "int")
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

# DDI
#human_microbiome_ddi = get_item(config_items, "DDI", "human_microbiome_ddi", "string")

# MSP
tshld_unclassified = 0.10	# the minimum percentile of unclassified in one MSP
tshld_diff = 0.50			# the minimum difference between most and second dominant taxa in the MSP
tshld_lca = 0.80			# the minimum consistency for LCA calculattion
taxa_final = "Rep"			# the source of taxa for one protein family, representatives vs. LCA
mspminer = os.path.join(config_directory, "MSPminer_setting.cfg")

# abundance
normalize = get_item(config_items, "abundance", "normalize", "string") 		# the method for abundance normalization
abundance_detection_level = get_item(config_items, "abundance", "abundance_detection_level", "float")

# maaslin2
maaslin2_dir =  os.path.join(abundance_dir, "DA", "maaslin2_output/") 
phenotype = get_item(config_items, "maaslin2", "phenotype", "string")
phenotype = re.sub("\"", "", phenotype)
phenotype = phenotype.split(";")
contrast_status = get_item(config_items, "maaslin2", "case_control_status", "string")
contrast_status = re.sub("\"", "", contrast_status)
tmp = contrast_status.split(";")
contrast_status = {}
for item in tmp:
	tmp1 = item.split(":")
	if len(tmp1) > 1:
		contrast_status[tmp1[0]] = tmp1[1]
ref_status = get_item(config_items, "maaslin2", "flag_ref", "string")
ref_status = re.sub("\"", "", ref_status)
tmp = ref_status.split(";")
ref_status = {}
for item in tmp:
	tmp1 = item.split(":")
	if len(tmp1) > 1:
		ref_status[tmp1[0]] = tmp1[1]
tshld_prevalence = get_item(config_items, "maaslin2", "tshld_prevalence", "float")
tshld_qvalue = get_item(config_items, "maaslin2", "tshld_qvalue", "float")
effect_size = get_item(config_items, "maaslin2", "effect_size", "string")
if effect_size == "log(fc)":
	effect_size = "log(FC)"
nested_effects = "none"
maaslin2_cmmd = get_item(config_items, "maaslin2", "maaslin2_cmmd", "string")
maaslin2_utils = os.path.join(metawibele_install_directory, "Rscripts", "maaslin2_utils.R")
pcl_utils = os.path.join(metawibele_install_directory, "Rscripts", "pcl_utils.R")
transpose_cmmd = "metawibele_transpose" 
min_abundance = get_item(config_items, "maaslin2", "min_abundance", "float")
min_prevalence = get_item(config_items, "maaslin2", "min_prevalence", "float")
#min_variance = get_item(config_items, "maaslin2", "min_variance", "float")
max_significance = get_item(config_items, "maaslin2", "max_significance", "float")
normalization = get_item(config_items, "maaslin2", "normalization", "string")
transform = get_item(config_items, "maaslin2", "transform", "string")
analysis_method = get_item(config_items, "maaslin2", "analysis_method", "string")
fixed_effects = get_item(config_items, "maaslin2", "fixed_effects", "string")
random_effects = get_item(config_items, "maaslin2", "random_effects", "string")
correction = get_item(config_items, "maaslin2", "correction", "string")
standardize = get_item(config_items, "maaslin2", "standardize", "string")
plot_heatmap = get_item(config_items, "maaslin2", "plot_heatmap", "string")
heatmap_first_n = get_item(config_items, "maaslin2", "heatmap_first_n", "string")
plot_scatter = get_item(config_items, "maaslin2", "plot_scatter", "string")
maaslin2_cores = get_item(config_items, "maaslin2", "maaslin2_cores", "int")
maaslin2_cmmd_opts = ["--min_abundance", min_abundance, "--min_prevalence", min_prevalence, "--max_significance", max_significance, "--normalization", normalization,  "--transform", transform, "--analysis_method", analysis_method, "--cores", maaslin2_cores, "--fixed_effects", fixed_effects, "--random_effects", random_effects, "--correction", correction, "--standardize", standardize, "--plot_heatmap", plot_heatmap, "--heatmap_first_n", heatmap_first_n, "--plot_scatter", plot_scatter]


## assembly workflow
# bowtie2 options and threshold
bowtie2_large_index_threshold = 4000000000
bowtie2_index_ext_list = [".1.bt2",".2.bt2",".3.bt2",".4.bt2",
    ".rev.1.bt2",".rev.2.bt2"]
bowtie2_large_index_ext = ".1.bt2l"
bowtie2_version = {
    "flag" : "--version",
    "major" : 2,
    "minor" : 2,
    "line" : 0,
    "column" : 2}
bowtie2_align_opts = ["--threads", threads] # "--threads "+str(threads)
bowtie2_index_name = "_bowtie2_index"


if __name__=='__main__':
	pass

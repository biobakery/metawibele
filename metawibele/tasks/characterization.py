"""
MetaWIBELE Workflow: characterization module
A collection of tasks for characterization workflow

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
import subprocess
import itertools
import re

from anadama2.tracked import TrackedExecutable, TrackedDirectory

# import the utilities functions and config settings from MetaWIBELE
try:
	from metawibele import utilities, config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the MetaWIBELE python package." +
		         " Please check your install.")


def clustering (workflow, gene_catalog_seq, threads, output_folder, protein_family, protein_family_seq):
	"""
	This set of tasks will run clustering on the input files provided.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		gene_catalog_seq (fasta file): Fasta file of amino acid sequences for gene catalogs.
		threads (int): The number of threads/cores for clustering to use.
		output_folder (string): The path of the output folder.
		protein_family (protein family file): This extending fasta file for clustering information
		protein_family_seq (protein family file): This fasta file of amino acid sequences for protein families.

	Requires:
		CD-hit version 4.7: A program for clustering and comparing protein or nucleotide sequences
		gene catalog sequences

	Returns:
		string: the name of protein families file.

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add clustering tasks
		protein_family, protein_family_output_folder = characterization.clustering(workflow, gene_catalog_seq,
                                                                           args.threads,
                                                                           output_dir,
                                                                           protein_family,
                                                                           protein_family_seq)
		# run the workflow
		workflow.go()
	"""

	print("clustering")
	
	# get the clustering output files
	main_folder = os.path.join(output_folder, "clustering")
	os.system("mkdir -p " + main_folder)
	tmps_dir = os.path.join (os.getcwd(), "temp")
	tmps_dir = os.path.abspath(tmps_dir)
	if not os.path.isdir(tmps_dir):
		os.system("mkdir -p " + tmps_dir)

	myname = os.path.basename(protein_family)
	clustering_output_cluster = os.path.join(main_folder, myname)
	myname = os.path.basename(protein_family_seq)
	clustering_output_seq = os.path.join(main_folder, myname)
	clustering_output_logs =  clustering_output_cluster + ".log"
	time_equation = config.time  # xxx hours defined in global config
	mem_equation = config.memory  # xxx GB defined in global config

	# get clustering option
	if config.cd_hit_prot_opts is None:
		optional_arguments = ""
	else:
		optional_arguments = ""
		tmp = config.cd_hit_prot_opts.split()
		for i in tmp:
			i = re.sub("\"", "", i)
			optional_arguments = optional_arguments + " " + str(i)
	if not threads is None:
		optional_arguments = "-T " + str(threads) + " " + optional_arguments

	# create tasks to run clustering
	myraw_cluster = clustering_output_seq + ".clstr"
	workflow.add_task(
			"cd-hit -i [depends[0]] " + optional_arguments + " -o [targets[0]] > [args[0]] 2>&1",
			depends = [gene_catalog_seq, TrackedExecutable("cd-hit")],
			targets = [clustering_output_seq, myraw_cluster],
			args = [clustering_output_logs],
			cores = threads,
			name = "clustering-proteins")

	workflow.add_task("ln -fs [depends[0]] [targets[0]]",
	        depends = [clustering_output_seq],
	        targets = [protein_family_seq],
	        cores = 1,
	        name = "clustering-proteins")
	workflow.add_task("ln -fs [depends[0]] [targets[0]]",
			depends = [clustering_output_seq],
	        targets = [os.path.join(tmps_dir, os.path.basename(protein_family_seq))],
	        cores = 1,
	        name = "ln__clustering-proteins")

	workflow.add_task(
			"metawibele_extract_cluster -c [depends[0]] -o [targets[0]] >> [args[0]] 2>&1",
			depends = [myraw_cluster, TrackedExecutable("metawibele_extract_cluster")],
			targets = [clustering_output_cluster],
			args = [clustering_output_logs],
			cores = 1,
			name = "extract_cluster_CD-hit")

	workflow.add_task("ln -fs [depends[0]] [targets[0]]",
			depends = [clustering_output_cluster],
	        targets = [protein_family],
	        cores = 1,
	        name = "ln__clustering-proteins")
	
	workflow.add_task("ln -fs [depends[0]] [targets[0]]",
			depends = [clustering_output_cluster],
	        targets = [os.path.join(tmps_dir, os.path.basename(protein_family))],
	        cores = 1,
	        name = "ln__clustering-proteins")

	return clustering_output_cluster, main_folder


def global_homology_annotation (workflow, family_conf, gene_catalog_seq,
                               threads, output_folder, uniref_taxonomy_family, uniref_taxonomy,
                               protein_family_ann_list, protein_ann_list, protein_family, protein_family_seq):
	"""
	This set of tasks will run annotations based on homologies to known proteins.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		family_conf: Methods configuration for homology-based annotation
		gene_catalog_seq (fasta file): Fasta file of amino acid sequences for gene catalogs.
		protein_family (protein family file): This extending fasta file for clustering information
		threads (int): The number of threads/cores for clustering to use.
		output_folder (string): The path of the output folder.
		uniref_annotation: The taxonomy annotation based on uniref90 database.

	Requires:
		diamond version 0.9.5:  a high-throughput program for aligning DNA reads or protein sequences against a protein reference database
		the amino acid sequences for gene catalogs
		the clustering information for protein families

	Returns:
		list: list of annotation files for protein/protein families directly extracted from UniRef90 database

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add global_homology_annotation tasks
		myprotein_family_ann, myprotein_ann, homology_output_folder = characterization.global_homology_annotation (workflow, family_conf, gene_catalog_seq,
	                                                                                                          protein_family,
	                                                                                                          args.threads,
	                                                                                                          output_dir, uniref_annotation)
		# run the workflow
		workflow.go()
	"""
	
	print("global_homology_annotation")

	# define the annotation output files
	main_folder = os.path.join(output_folder, "global_homology_annotation")
	os.system("mkdir -p " + main_folder)
	tmps_dir = os.path.join (os.getcwd(), "temp")
	tmps_dir = os.path.abspath(tmps_dir)
	if not os.path.isdir(tmps_dir):
		os.system("mkdir -p " + tmps_dir)

	myname = os.path.basename(gene_catalog_seq)
	myname = re.sub(".tsv", "", myname)
	annotation_stat = os.path.join(main_folder, myname + ".uniref90.stat.tsv")
	uniref_ann_protein = os.path.join(main_folder, config.basename + "_UniRef90_protein.tsv")
	uniref_ann_family = os.path.join(main_folder,config.basename + "_UniRef90_proteinfamilies.detail.tsv")
	uniref_ann = os.path.join(main_folder, config.basename + "_UniRef90_proteinfamilies.ORF.detail.tsv")
	myname = os.path.basename(uniref_taxonomy_family)
	uniref_taxa_family = os.path.join(main_folder, myname)
	uniref_taxa = re.sub("proteinfamilies", "protein", uniref_taxa_family)
	uniref_log0 = os.path.join(main_folder, "uniref90.mapping.log")
	uniref_log1 =  os.path.join(main_folder, "uniref90.stat.log")
	uniref_log2 = os.path.join(main_folder, "uniref90_protein.log")
	uniref_log3 = os.path.join(main_folder, "uniref90_proteinfamilies.log")
	uniref_log4 = os.path.join(main_folder, "uniref90_annotation.log")
	uniref_log5 = os.path.join(main_folder, "uniref90_family_annotation.log")
	myprotein_family_ann = {}
	myprotein_ann = {}
	myprotein_ann_family = {}
	time_equation = config.time  # xxx hours defined in global config
	mem_equation = config.memory  # xxx GB defined in global config

	# mapping to UniRef database
	if family_conf["uniref"].lower() == "yes":
		myname = os.path.basename(gene_catalog_seq)
		myhit = re.sub(".uniref90.stat.tsv", ".uniref90.hits", annotation_stat)
		link_cmd = "ln -fs " + gene_catalog_seq + " " + main_folder + "/" + myname
		os.system(link_cmd)
		workflow.add_task_gridable(
				link_cmd + " && metawibele_uniref_annotator [depends[0]] --seqtype prot --uniref90db [depends[1]] --diamond-options \"--threads [args[0]]\" --temp [args[1]] >[args[2]] 2>&1",
				depends = [os.path.join(main_folder, myname), config.uniref_dmnd, TrackedExecutable("metawibele_uniref_annotator")],
				targets = [myhit],
				args = [threads, main_folder, uniref_log0],
				cores = threads,
				time = time_equation,
				mem = mem_equation,
				name = "uniref_annotator")

		# extract mapping info
		workflow.add_task(
				"metawibele_uniref_annotator_stat -s [depends[0]] -d [depends[1]] "
				"-o [targets[0]] >[args[0]] 2>&1",
				depends = [gene_catalog_seq, myhit, TrackedExecutable("metawibele_uniref_annotator_stat")],
				targets = [annotation_stat],
				args = [uniref_log1],
				cores = 1,
				name = "uniref_annotator_stat")

		# uniRef annotation for each ORF
		#protein_family_tmp = os.path.join(tmps_dir, os.path.basename(protein_family))
		#protein_family_seq_tmp = os.path.join(tmps_dir, os.path.basename(protein_family_seq))
		protein_family_tmp = protein_family
		protein_family_seq_tmp = protein_family_seq
		workflow.add_task(
				"metawibele_uniref_protein -m [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [annotation_stat, protein_family_tmp, protein_family_seq_tmp, TrackedExecutable("metawibele_uniref_protein")],
				targets = [uniref_ann_protein],
				args = [uniref_log2],
				cores = 1,
				name = "uniref_protein")

		# uniRef annotation for each protein family
		workflow.add_task(
				"metawibele_uniref_protein_family -u [depends[0]] -m [depends[1]] -f centroid -o [targets[0]] >[args[0]] 2>&1",
				depends = [uniref_ann_protein, annotation_stat, protein_family_tmp, protein_family_seq_tmp, TrackedExecutable("metawibele_uniref_protein_family")],
				targets = [uniref_ann_family, uniref_ann],
				args = [uniref_log3],
				cores = 1,
				name = "uniref_protein_family")


		# uniref annotation for homology and taxonomy
		workflow.add_task(
				"metawibele_summary_protein_uniref_annotation -a [depends[0]] -m [depends[1]] -t Rep -o [targets[0]] >[args[0]] 2>&1",
				depends = [uniref_ann, annotation_stat, TrackedExecutable("metawibele_summary_protein_uniref_annotation")],
				targets = [uniref_taxa],
				args = [uniref_log4],
				cores = 1,
				name = "summary_protein_uniref_annotation")

		workflow.add_task(
				"metawibele_summary_protein_family_uniref_annotation -a [depends[0]] -m [depends[1]] -t Rep -o [targets[0]] >[args[0]] 2>&1",
				depends = [uniref_ann_family, annotation_stat, TrackedExecutable("metawibele_summary_protein_family_uniref_annotation")],
				targets = [uniref_taxa_family],
				args = [uniref_log5],
				cores = 1,
				name = "summary_protein_family_uniref_annotation")

		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
		        depends = [uniref_taxa_family],
	            targets = [uniref_taxonomy_family],
		        cores = 1,
	            name = "ln__uniref_taxa_family")
		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
	            depends = [uniref_taxa],
	            targets = [uniref_taxonomy],
		        cores = 1,
	            name = "ln__uniref_taxa")
		
		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
		        depends = [uniref_taxa_family],
	            targets = [os.path.join(tmps_dir, os.path.basename(uniref_taxonomy_family))],
		        cores = 1,
	            name = "ln__uniref_taxa_family")
		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
	            depends = [uniref_taxa],
	            targets = [os.path.join(tmps_dir, os.path.basename(uniref_taxonomy))],
		        cores = 1,
	            name = "ln__uniref_taxa")
		
		myprotein_family_ann[uniref_ann_family] = ""
		myprotein_ann[uniref_ann] = ""
	# if uniref annotation

	if len(myprotein_family_ann.keys()) > 0:
		for myfile in myprotein_family_ann.keys():
			protein_family_ann_list[myfile] = ""
	if len(myprotein_ann.keys()) > 0:
		for myfile in myprotein_ann.keys():
			protein_ann_list[myfile] = ""

	return myprotein_ann_family, myprotein_ann, main_folder


def domain_motif_annotation (workflow, domain_motif_conf, gene_catalog_seq,
                               split_number, threads, output_folder,
                               protein_family_ann_list, protein_ann_list, protein_family, protein_family_seq):
	"""
	This set of tasks will run annotations predicted by sequences.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		domain_motif_conf: Methods configuration for sequence-based annotation
		gene_catalog_seq (fasta file): Fasta file of amino acid sequences for gene catalogs.
		protein_family (protein family file): This extending fasta file for clustering information
		threads (int): The number of threads/cores for clustering to use.
		output_folder (string): The path of the output folder.

	Requires:
		InterProScan version 5.31-70.0: is the software package that allows sequences (protein and nucleic) to be scanned against InterPro's signatures. Signatures are predictive models, provided by several different databases
		PSORTb version 3.0: protein subcellular localization prediction
		signalp version 4.1: predicts the presence of signal peptides and the location of their cleavage sites in proteins
		TMHMM version 2.0c: Prediction of transmembrane helices in proteins
		Phobius version 1.01: A combined transmembrane topology and signal peptide predictor
		the amino acid sequences for gene catalogs
		the clustering information for protein families

	Returns:
		list: list of annotation files for protein/protein families predicted based on sequences

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add domain_motif_annotation tasks
		uniref_annotation, protein_family_ann, protein_ann, sequence_output_folder = characterization.domain_motif_annotation (workflow, domain_motif_conf, gene_catalog_seq,
	                                                                                                                         protein_family,
	                                                                                                                         args.threads,
	                                                                                                                         output_dir)

		# run the workflow
		workflow.go()
	"""
	
	print("domain_motif_annotation")

	# define the annotation output files
	main_folder = os.path.join(output_folder, "domain_motif_annotation")
	interpro = os.path.join(main_folder, "InterProScan")
	if len(interpro) >200:
		# debug
		sys.stderr.write("The path for InterProScan output has more than 230 characters it throws the error you are getting! Change your path to home folder.")
		interpro = "~/"
	psortb = os.path.join(main_folder, "PSORTb")
	tmps_dir = os.path.join (os.getcwd(), "temp")
	tmps_dir = os.path.abspath(tmps_dir)
	if not os.path.isdir(tmps_dir):
		os.system("mkdir -p " + tmps_dir)

	commands = "mkdir -p " + main_folder
	commands = commands + "; " + "mkdir -p " + interpro
	commands = commands + "; " + "mkdir -p " + psortb
	os.system(commands)

	myprefix = os.path.basename(gene_catalog_seq)
	myprefix = re.sub(".tsv", "", myprefix)
	signalp_ann_family = os.path.join(main_folder, config.basename + "_SignalP_proteinfamilies.detail.tsv")
	signalp_ann = os.path.join(main_folder, config.basename + "_SignalP_proteinfamilies.ORF.detail.tsv")
	tmhmm_ann_family = os.path.join(main_folder, config.basename + "_TMHMM_proteinfamilies.detail.tsv")
	tmhmm_ann = os.path.join(main_folder, config.basename + "_TMHMM_proteinfamilies.ORF.detail.tsv")
	phobius_ann_family = os.path.join(main_folder, config.basename + "_Phobius_proteinfamilies.detail.tsv")
	phobius_ann = os.path.join(main_folder, config.basename + "_Phobius_proteinfamilies.ORF.detail.tsv")
	denovo_ann_family = os.path.join(main_folder, config.basename + "_Denovo_proteinfamilies.detail.tsv")
	denovo_ann = os.path.join(main_folder, config.basename + "_Denovo_proteinfamilies.ORF.detail.tsv")
	denovo_signal_ann_family = os.path.join(main_folder, config.basename + "_Denovo_proteinfamilies.signaling.detail.tsv")
	denovo_signal_ann = os.path.join(main_folder, config.basename + "_Denovo_proteinfamilies.ORF.signaling.detail.tsv")
	denovo_trans_ann_family = os.path.join(main_folder, config.basename + "_Denovo_proteinfamilies.transmembrane.detail.tsv")
	denovo_trans_ann = os.path.join(main_folder, config.basename + "_Denovo_proteinfamilies.ORF.transmembrane.detail.tsv")
	pfam_ann_family = os.path.join(main_folder, config.basename + "_Pfam_proteinfamilies.detail.tsv")
	pfam_ann = os.path.join(main_folder, config.basename + "_Pfam_proteinfamilies.ORF.detail.tsv")
	domine_ann_family = os.path.join(main_folder, config.basename + "_DOMINE_proteinfamilies.detail.tsv")
	domine_ann = os.path.join(main_folder, config.basename + "_DOMINE_proteinfamilies.ORF.detail.tsv")
	domine_ann_family_all = os.path.join(main_folder, config.basename + "_DOMINE_proteinfamilies.all.detail.tsv")
	domine_ann_all = os.path.join(main_folder, config.basename + "_DOMINE_proteinfamilies.ORF.all.detail.tsv")
	SIFTS_ann_family = os.path.join(main_folder, config.basename + "_SIFTS_proteinfamilies.detail.tsv")
	SIFTS_ann = os.path.join(main_folder, config.basename + "_SIFTS_proteinfamilies.ORF.detail.tsv")
	ExpAtlas_ann_family = os.path.join(main_folder, config.basename + "_ExpAtlas_proteinfamilies.detail.tsv")
	ExpAtlas_ann = os.path.join(main_folder, config.basename + "_ExpAtlas_proteinfamilies.ORF.detail.tsv")
	pfam2go_ann_family = os.path.join(main_folder, config.basename + "_Pfam2GO_proteinfamilies.detail.tsv")
	pfam2go_ann = os.path.join(main_folder, config.basename + "_Pfam2GO_proteinfamilies.ORF.detail.tsv")
	interProScan_ann_family = os.path.join(main_folder, config.basename + "_InterProScan_proteinfamilies.detail.tsv")
	interProScan_ann = os.path.join(main_folder, config.basename + "_InterProScan_proteinfamilies.ORF.detail.tsv")
	psortb_ann_family = os.path.join(main_folder, config.basename + "_PSORTb_proteinfamilies.detail.tsv")
	psortb_ann = os.path.join(main_folder, config.basename + "_PSORTb_proteinfamilies.ORF.detail.tsv")
	myprotein_family_ann = {}
	myprotein_ann = {}
	time_equation = config.time  # xxx hours defined in global config
	mem_equation = config.memory  # xxx GB defined in global config
	
	# split files as inputs
	file_list_file = os.path.join(tmps_dir, "split_files.list")
	split_list_file = os.path.join(tmps_dir, "split.list")
	utilities.split_fasta_file(gene_catalog_seq, split_number, myprefix, tmps_dir, file_list_file, split_list_file)
	checking = 1
	while checking:
		if utilities.is_file_exist(file_list_file):
			file_list = []
			open_file = open(file_list_file, "r")
			for line in open_file:
				line = line.strip()
				file_list.append(line)
			open_file.close()
			checking = 0
	checking = 1
	while checking:
		if utilities.is_file_exist(split_list_file):
			split_list = []
			open_file = open(split_list_file, "r")
			for line in open_file:
				line = line.strip()
				split_list.append(line)
			open_file.close()
			checking = 0

	interpro_list1 = []
	interpro_list = []
	pfam_list = []
	for myfile, mysplit in zip(file_list, split_list):
		myout_dir = os.path.join(interpro, mysplit)
		myout = os.path.join(myout_dir, mysplit + ".interproscan.txt")
		interpro_list1.append(myout)
		myfile_new = os.path.join(interpro, os.path.basename(myout))
		interpro_list.append(myfile_new)
		myfile1 = re.sub("interproscan.txt", "interpro.PfamDomain.tsv", myfile_new)
		pfam_list.append(myfile1)
	
	## run InterProScan per each protein
	if domain_motif_conf["interproscan"].lower() == "yes":
		for myfile, mysplit in zip(file_list, split_list):
			myout_dir = os.path.join(interpro, mysplit)
			os.system('mkdir -p ' + myout_dir)
			myout = os.path.join(myout_dir, mysplit + ".interproscan.txt")
			mylog = os.path.join(myout_dir, os.path.basename(myfile))
			mylog = re.sub(".fasta", ".log", mylog)
			myerr = re.sub(".log", ".err", mylog)

			workflow.add_task_gridable(
					"mkdir -p [args[2]] && metawibele_interproscan_annotator --split-file [args[0]] --threads [args[1]] -o [args[2]] -i [depends[0]] > [args[3]] 2> [args[4]] ",
					depends = [myfile, TrackedExecutable("metawibele_interproscan_annotator")],
					targets = [myout],
					args = [mysplit, threads, myout_dir, mylog, myerr],
					cores = threads,   
					time = time_equation,
					mem = mem_equation,
					name = utilities.name_task(mysplit, "interproscan"))
		
		for myfile in interpro_list1:
			myname = os.path.basename(myfile)
			myfile_new = os.path.join(interpro,  myname)
			workflow.add_task(
					"ln -fs [depends[0]] [targets[0]]",
					depends = [myfile],
					targets = [myfile_new],
					cores = 1,
					name = utilities.name_task(myname, "ln"))

		## InterProScan annotation for each ORF
		mylog = os.path.join(interpro, "interproscan.extract.log")
		myout10 = []
		for myfile in interpro_list:
			myfile1 = re.sub("interproscan.txt", "signalp.signaling.tsv", myfile)
			myout10.append(myfile1)
			myfile1 = re.sub("interproscan.txt", "tmhmm.transmembrane.tsv", myfile)
			myout10.append(myfile1)
			myfile1 = re.sub("interproscan.txt", "phobius.signaling.tsv", myfile)
			myout10.append(myfile1)
			myfile1 = re.sub("interproscan.txt", "phobius.transmembrane.tsv", myfile)
			myout10.append(myfile1)
			myfile1 = re.sub("interproscan.txt", "interpro.PfamDomain.tsv", myfile)
			myout10.append(myfile1)
		workflow.add_task(
				"metawibele_interproscan_protein -e [args[0]] -p [args[1]] >[args[2]] 2>&1",
				depends = utilities.add_to_list(interpro_list, TrackedExecutable("metawibele_interproscan_protein")),
				targets = myout10,
				args = ["interproscan.txt", interpro, mylog],
				cores = 1,
				name = "interproscan_protein")

		## InterProScan annotation for each family
		mylog = re.sub(".tsv", ".log", signalp_ann_family)
		myout20 = myout10
		myout20.append(protein_family)
		myout20.append(protein_family_seq)
		workflow.add_task(
				"metawibele_interproscan_signalp_protein_family -e [args[0]] -p [args[1]] -a consistency -o [targets[0]] >[args[2]] 2>&1",
				depends = utilities.add_to_list(myout20, TrackedExecutable("metawibele_interproscan_signalp_protein_family")), 
				targets = [signalp_ann_family, signalp_ann],
				args = ["signalp.signaling.tsv", interpro, mylog],
				cores = 1,
				name = "interproscan_signalp_protein_family")

		mylog = re.sub(".tsv", ".log", tmhmm_ann_family)
		workflow.add_task(
				"metawibele_interproscan_tmhmm_protein_family -e [args[0]] -p [args[1]] -a consistency -o [targets[0]] >[args[2]] 2>&1",
				depends = utilities.add_to_list(myout20, TrackedExecutable("metawibele_interproscan_tmhmm_protein_family")),
				targets = [tmhmm_ann_family, tmhmm_ann],
				args = ["tmhmm.transmembrane.tsv", interpro, mylog],
				cores = 1,
				name = "interproscan_tmhmm_protein_family")
		
		myfile1 = re.sub(".detail.tsv", ".signaling.detail.tsv", phobius_ann_family)
		myfile2 = re.sub(".detail.tsv", ".transmembrane.detail.tsv", phobius_ann_family)
		myfile3 = re.sub(".detail.tsv", ".signaling.detail.tsv", phobius_ann)
		myfile4 = re.sub(".detail.tsv", ".transmembrane.detail.tsv", phobius_ann)
		mylog = re.sub(".tsv", ".log", phobius_ann_family)
		workflow.add_task(
				"metawibele_interproscan_phobius_protein_family -e [args[1]] -p [args[2]] -a consistency -o [args[0]] >[args[3]] 2>&1",
				depends = utilities.add_to_list(myout20, TrackedExecutable("metawibele_interproscan_phobius_protein_family")),
				targets = [myfile1, myfile2, myfile3, myfile4],
				args = [phobius_ann_family, "phobius.signaling.tsv", interpro, mylog],
				cores = 1,
				name = "interproscan_phobius_protein_family")

		mylog = re.sub(".tsv", ".log", pfam_ann_family)
		workflow.add_task(
				"metawibele_interproscan_pfam_protein_family -e [args[0]] -p [args[1]] -a consistency -o [targets[0]] >[args[2]] 2>&1",
				depends = utilities.add_to_list(myout20, TrackedExecutable("metawibele_interproscan_pfam_protein_family")),
				targets = [pfam_ann_family, pfam_ann],
				args = ["interpro.PfamDomain.tsv", interpro, mylog],
				cores = 1,
				name = "interproscan_pfam_protein_family")

		mylog = re.sub(".tsv", ".log", interProScan_ann_family)
		workflow.add_task(
				"metawibele_interproscan_protein_family -e [args[0]] -p [args[1]] -a consistency -o [targets[0]] >[args[2]] 2>&1",
				depends = utilities.add_to_list(myout20, TrackedExecutable("metawibele_interproscan_protein_family")),
				targets = [interProScan_ann_family, interProScan_ann],
				args = ["interproscan.txt", interpro, mylog],
				cores = 1,
				name = "interproscan_protein_family")
		
		## build de novo sigaling/transmembrane set
		mylog = re.sub(".tsv", ".log", denovo_ann_family)
		myfile1 = re.sub(".detail.tsv", ".signaling.detail.tsv", phobius_ann_family)
		myfile2 = re.sub(".detail.tsv", ".transmembrane.detail.tsv", phobius_ann_family)
		workflow.add_task(
				"metawibele_denovo_TM_SP -a [depends[0]] -b [depends[1]] -c [depends[2]] -d [depends[3]] -o [args[0]] >[args[1]] 2>&1",
				depends = [myfile1, myfile2, signalp_ann_family, tmhmm_ann_family, protein_family, protein_family_seq, TrackedExecutable("metawibele_denovo_TM_SP")],
				targets = [denovo_signal_ann_family, denovo_trans_ann_family],
				args = [denovo_ann_family, mylog],
				cores = 1,
				name = "denovo_TM_SP")

		mylog = re.sub(".tsv", ".log", denovo_ann)
		myfile1 = re.sub(".detail.tsv", ".signaling.detail.tsv", phobius_ann)
		myfile2 = re.sub(".detail.tsv", ".transmembrane.detail.tsv", phobius_ann)
		workflow.add_task(
				"metawibele_denovo_TM_SP -a [depends[0]] -b [depends[1]] -c [depends[2]] -d [depends[3]] -o [args[0]] > [args[1]] 2>&1",
				depends = [myfile1, myfile2, signalp_ann, tmhmm_ann, TrackedExecutable("metawibele_denovo_TM_SP")],
				targets = [denovo_signal_ann, denovo_trans_ann],
				args = [denovo_ann, mylog],
				cores = 1,
				name = "denovo_TM_SP")

		myprotein_family_ann[interProScan_ann_family] = ""
		myprotein_ann[interProScan_ann] = ""
		myprotein_ann[denovo_signal_ann] = ""
		myprotein_ann[denovo_trans_ann] = ""
		myprotein_family_ann[denovo_signal_ann_family] = ""
		myprotein_family_ann[denovo_trans_ann_family] = ""
	# if InterProScan annotation


	## Pfam2GO info
	if domain_motif_conf["pfam2go"].lower() == "yes":
		if domain_motif_conf["interproscan"].lower() == "no":
			sys.stderr.write("WARNING: Please make sure that domain annotations have been assigned by interproscan. Otherwise, errors might be caused!")
		mylog = re.sub(".tsv", ".log", pfam2go_ann_family)
		workflow.add_task(
				"metawibele_pfam2go -i [depends[0]] -o [targets[0]] > [args[0]] 2>&1",
				depends = [pfam_ann_family, protein_family, protein_family_seq, TrackedExecutable("metawibele_pfam2go")],
				targets = [pfam2go_ann_family],
				args = [mylog],
				cores = 1,
				name = "pfam2go")

		mylog = re.sub(".tsv", ".log", pfam2go_ann)
		workflow.add_task(
				"metawibele_pfam2go -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1 ",
				depends = [pfam_ann, TrackedExecutable("metawibele_pfam2go")],
				targets = [pfam2go_ann],
				args = [mylog],
				cores = 1,
				name = "pfam2go")

		myprotein_family_ann[pfam2go_ann_family] = ""
		myprotein_ann[pfam2go_ann] = ""
	# if Pfam2GO

	## DOMINE (domain-domain interaction)
	if domain_motif_conf["domine"].lower() == "yes":
		if domain_motif_conf["interproscan"].lower() == "no":
			sys.stderr.write("WARNING: Please make sure that domain annotations have been assigned by interproscan. Otherwise, errors might be caused!")
		# DDI annotation
		mylog = re.sub(".tsv", ".log", domine_ann)
		myout = []
		myout_all = []
		for myfile in interpro_list:
			myfile = re.sub("interproscan.txt", "interpro.DDI.tsv", myfile)
			myout.append(myfile)
			myfile = re.sub("interproscan.txt", "interpro.all.DDI.tsv", myfile)
			myout_all.append(myfile)
		
		workflow.add_task(
				"metawibele_ddi_DOMINE_protein -e [args[0]] -p [args[2]] -f [args[3]] -s [args[1]] >[args[4]] 2>&1 ",
				depends = utilities.add_to_list(pfam_list, TrackedExecutable("metawibele_ddi_DOMINE_protein")),
				targets = myout,
				args = ["interpro.PfamDomain.tsv", "interpro.DDI.tsv", interpro, "yes", mylog],
				cores = 1,
				name = "ddi_DOMINE_protein")
		
		workflow.add_task(
				"metawibele_ddi_DOMINE_protein -e [args[0]] -p [args[2]] -f [args[3]] -s [args[1]] >> [args[4]] 2>&1 ",
				depends = utilities.add_to_list(pfam_list, TrackedExecutable("metawibele_ddi_DOMINE_protein")),
				targets = myout_all,
				args = ["interpro.PfamDomain.tsv", "interpro.all.DDI.tsv", interpro, "no", mylog],
				cores = 1,
				name = "ddi_DOMINE_protein")

		mylog = re.sub(".tsv", ".log", domine_ann_family)
		domine_ann_family_raw = re.sub(".detail", "", domine_ann_family)
		domine_ann_raw = re.sub(".detail", "", domine_ann)
		myout2 = myout
		myout2.append(protein_family)
		myout2.append(protein_family_seq)
		myout3 = myout_all
		myout3.append(protein_family)
		myout3.append(protein_family_seq)
		workflow.add_task(
				"metawibele_ddi_DOMINE_protein_family -e [args[0]] -p [args[1]] -a consistency -o [targets[0]] >[args[2]] 2>&1 ",
				depends = utilities.add_to_list(myout2, TrackedExecutable("metawibele_ddi_DOMINE_protein_family")),
				targets = [domine_ann_family_raw, domine_ann_raw, domine_ann_family, domine_ann],
				args = ["interpro.DDI.tsv", interpro, mylog],
				cores = 1,
				name = "ddi_DOMINE_protein_family")
		
		mylog = re.sub(".tsv", ".log", domine_ann_family_all)
		domine_ann_family_all_raw = re.sub(".detail", "", domine_ann_family_all)
		domine_ann_all_raw = re.sub(".detail", "", domine_ann_all)
		workflow.add_task(
				"metawibele_ddi_DOMINE_protein_family -e [args[0]] -p [args[1]] -a consistency -o [targets[0]] >[args[2]] 2>&1 ",
				depends = utilities.add_to_list(myout3, TrackedExecutable("metawibele_ddi_DOMINE_protein_family")),
				targets = [domine_ann_family_all_raw, domine_ann_all_raw, domine_ann_family_all, domine_ann_all],
				args = ["interpro.all.DDI.tsv", interpro, mylog],
				cores = 1,
				name = "ddi_DOMINE_protein_family")

		myout1_2 = re.sub(".tsv", ".ann.tsv", domine_ann_family)
		mylog = re.sub(".tsv", ".log", myout1_2)
		workflow.add_task(
				"metawibele_ddi_DOMINE_ann -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [domine_ann_family_raw, protein_family, protein_family_seq, TrackedExecutable("metawibele_ddi_DOMINE_ann")],
				targets = [myout1_2],
				args = [mylog],
				cores = 1,
				name = "ddi_DOMINE_ann")

		myout2_2 = re.sub(".tsv", ".ann.tsv", domine_ann)
		mylog = re.sub(".tsv", ".log", myout2_2)
		workflow.add_task(
				"metawibele_ddi_DOMINE_ann -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [domine_ann_raw, TrackedExecutable("metawibele_ddi_DOMINE_ann")],
				targets = [myout2_2],
				args = [mylog],
				cores = 1,
				name = "ddi_DOMINE_ann")

		myprotein_family_ann[domine_ann_family] = ""
		myprotein_ann[domine_ann] = ""
		myprotein_family_ann[domine_ann_family_all] = ""
		myprotein_ann[domine_ann_all] = ""
	# if DDI

	# DDI + SIFTS annotation
	if domain_motif_conf["sifts"].lower() == "yes":
		if domain_motif_conf["interproscan"].lower() == "no":
			sys.stderr.write("WARNING: Please make sure that domain annotations have been assigned by interproscan. Otherwise, errors might be caused!")
		if domain_motif_conf["domine"].lower() == "no":
			sys.stderr.write("WARNING: Please make sure that DDIs have been assigned based on DOMINE. Otherwise, errors might be caused!")
		mylog = re.sub(".tsv", ".log", SIFTS_ann_family)
		myout1_1 = re.sub(".tsv", ".ann.tsv", domine_ann_family)
		myout2_1 = re.sub(".tsv", ".ann.tsv", domine_ann)
		workflow.add_task(
				"metawibele_ddi_DOMINE_SIFTS -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [myout1_1, protein_family, protein_family_seq, TrackedExecutable("metawibele_ddi_DOMINE_SIFTS")],
				targets = [SIFTS_ann_family],
				args = [mylog],
				cores = 1,
				name = "ddi_DOMINE_SIFTS")

		mylog = re.sub(".tsv", ".log", SIFTS_ann)
		workflow.add_task(
				"metawibele_ddi_DOMINE_SIFTS -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [myout2_1, TrackedExecutable("metawibele_ddi_DOMINE_SIFTS")],
				targets = [SIFTS_ann],
				args = [mylog],
				cores = 1,
				name = "ddi_DOMINE_SIFTS")

		myprotein_family_ann[SIFTS_ann_family] = ""
		myprotein_ann[SIFTS_ann] = ""
	# if SIFTS

	# DDI + human expression
	if domain_motif_conf["expatlas"].lower() == "yes": 
		if domain_motif_conf["interproscan"].lower() == "no":
			sys.stderr.write("WARNING: Please make sure that domain annotations have been assigned by interproscan. Otherwise, errors might be caused!")
		if domain_motif_conf["domine"].lower() == "no":
			sys.stderr.write("WARNING: Please make sure that DDIs have been assigned based on DOMINE. Otherwise, errors might be caused!")
		mylog = re.sub(".tsv", ".log", ExpAtlas_ann_family)
		myout1_3 = re.sub(".tsv", ".ann.tsv", domine_ann_family)
		myout2_3 = re.sub(".tsv", ".ann.tsv", domine_ann)
		workflow.add_task(
				"metawibele_ddi_DOMINE_ExpAtlas -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [myout1_3, protein_family, protein_family_seq, TrackedExecutable("metawibele_ddi_DOMINE_ExpAtlas")],
				targets = [ExpAtlas_ann_family],
				args = [mylog],
				cores = 1,
				name = "ddi_DOMINE_ExpAtlas")

		mylog = re.sub(".tsv", ".log", ExpAtlas_ann)
		workflow.add_task(
				"metawibele_ddi_DOMINE_ExpAtlas -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [myout2_3, TrackedExecutable("metawibele_ddi_DOMINE_ExpAtlas")],
				targets = [ExpAtlas_ann],
				args = [mylog],
				cores = 1,
				name = "ddi_DOMINE_ExpAtlas")

		myprotein_family_ann[ExpAtlas_ann_family] = ""
		myprotein_ann[ExpAtlas_ann] = ""
	# if ExpAtlas

	## run PSORTb
	if domain_motif_conf["psortb"].lower() == "yes":
		psortb_list1 = []
		for myfile, mysplit in zip(file_list, split_list):
			myout_dir = os.path.join(psortb, mysplit)
			os.system('mkdir -p ' + myout_dir)
			myout1_4 = os.path.join(myout_dir, mysplit + ".psortb.gram_positive.out.txt")
			psortb_list1.append(myout1_4)
			myout2_4 = os.path.join(myout_dir, mysplit + ".psortb.gram_negative.out.txt") 
			psortb_list1.append(myout2_4)
			myout3_4 = os.path.join(myout_dir, mysplit + ".psortb.archaea.out.txt") 
			psortb_list1.append(myout3_4)
			mylog = os.path.join(myout_dir, os.path.basename(myfile))
			mylog = re.sub(".fasta", ".log", mylog)
			myerr = re.sub(".log", ".err", mylog)

			workflow.add_task_gridable(
					"mkdir -p [args[2]] && metawibele_psortb_annotator --split-file [args[0]] --threads [args[1]] -o [args[2]] -i [depends[0]] > [args[3]] 2> [args[4]] ",
					depends = [myfile, TrackedExecutable("metawibele_psortb_annotator")],
					targets = [myout1_4, myout2_4, myout3_4],
					args = [mysplit, threads, myout_dir, mylog, myerr],
					cores = threads,
					time = time_equation,
					mem = mem_equation,
					name = utilities.name_task(mysplit, "psortb"))
		# foreach split file 

		psortb_list = []
		myout5 = [] 
		for myfile in psortb_list1:
			myname = os.path.basename(myfile)
			myfile_new = os.path.join(psortb, myname)
			psortb_list.append(myfile_new) 
			myfile1 = re.sub("psortb.gram_positive.out.txt", "psortb.gram_positive.out.location.tsv", myfile_new)
			myfile1 = re.sub("psortb.gram_negative.out.txt", "psortb.gram_negative.out.location.tsv", myfile_new)
			myfile1 = re.sub("psortb.archaea.out.txt", "psortb.archaea.out.location.tsv", myfile_new)
			myout5.append(myfile1)
			workflow.add_task(
					"ln -fs [depends[0]] [targets[0]]",
					depends = [myfile],
					targets = [myfile_new],
					cores = 1,
					name = utilities.name_task(myname, "ln"))

		## PSORTb annotation for each ORF
		mylog = os.path.join(psortb, "psortb.extract.log")
		workflow.add_task(
				"metawibele_psortb_protein -e [args[0]] -p [args[1]] >[args[2]] 2>&1",
				depends = utilities.add_to_list(psortb_list, TrackedExecutable("metawibele_psortb_protein")),	
				targets = myout5,
				args = ["psortb.gram_positive.out.txt", psortb, mylog],
				cores = 1,
				name = "psortb_protein")

		## PSORTb annotation for each family
		mylog = re.sub(".tsv", ".log", psortb_ann_family)
		myout2_5 = myout5
		myout2_5.append(protein_family)
		myout2_5.append(protein_family_seq)
		workflow.add_task(
				"metawibele_psortb_protein_family -e [args[0]] -p [args[1]] -a consistency -o [targets[0]] >[args[2]] 2>&1",
				depends = utilities.add_to_list(myout2_5, TrackedExecutable("metawibele_psortb_protein_family")),
				targets = [psortb_ann_family, psortb_ann],
				args = ["psortb.gram_positive.out.location.tsv", psortb, mylog],
				cores = 1,
				name = "psortb_protein_family")

		myprotein_family_ann[psortb_ann_family] = ""
		myprotein_ann[psortb_ann] = ""
	# if PSORTb

	if len(myprotein_family_ann.keys()) > 0:
		for myfile in myprotein_family_ann.keys():
			protein_family_ann_list[myfile] = ""
	if len(myprotein_ann.keys()) > 0:
		for myfile in myprotein_ann.keys():
			protein_ann_list[myfile] = ""

	return myprotein_family_ann, myprotein_ann, main_folder


def abundance_annotation (workflow, abundance_conf, gene_catalog_seq, gene_catalog_count,
                                uniref_taxonomy_family, uniref_taxonomy, split_number,
                                threads, output_folder, protein_family, protein_family_relab, taxonomy_annotation_family, taxonomy_annotation,
                                protein_family_ann_list, protein_ann_list):
	"""
	This set of tasks will run annotations predicted by abundance.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		abundance_conf: Methods configuration for abundance-based annotation
		gene_catalog_seq (fasta file): Fasta file of amino acid sequences for gene catalogs.
		gene_catalog_count (count table): reads count for gene catalogs.
		protein_family (protein family file): This extending fasta file for clustering information
		threads (int): The number of threads/cores for clustering to use.
		output_folder (string): The path of the output folder.

	Requires:
		MSPminer V2: a tool to create Metagenomic Species Pan-genomes from genes abundance profiles
		MaAsLin2: a comprehensive R package for efficiently determining multivariable association between microbial meta'omic features and clinical metadata
		the counts table for gene catalogs
		the clustering information for protein families
		the uniref90 taxonomic annotations for protein families

	Returns:
		list: list of annotation files for protein/protein families predicted based on sequences

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add abundance_annotation tasks
		taxonomy_annotation, protein_family_ann, protein_ann, abundance_output_folder = characterization.abundance_annotation (workflow, abundance_conf, gene_catalog_count,
	                                                                                                                             protein_family,
	                                                                                                                             uniref_taxonomy,
	                                                                                                                             split_number,
	                                                                                                                             args.threads,
	                                                                                                                             output_dir,
	                                                                                                                             protein_family_relab, taxonomy_annotation)

		# run the workflow
		workflow.go()
	"""

	print("abundance_annotation")
	
	# define the annotation output files
	main_folder = os.path.join(output_folder, "abundance_annotation")
	os.system("mkdir -p " + main_folder)
	tmps_dir = os.path.join (os.getcwd(), "temp")
	tmps_dir = os.path.abspath(tmps_dir)
	if not os.path.isdir(tmps_dir):
		os.system("mkdir -p " + tmps_dir)

	myname = os.path.basename(gene_catalog_count)
	myprefix = re.sub(".tsv", "", myname)
	count_file= os.path.join(main_folder,  myprefix + ".refined.tsv")
	mymsp = os.path.join(main_folder, config.basename + "_MSPminer_msp.tsv")
	mymsp_uniref = os.path.join(main_folder, config.basename + "_MSPminer_msp.uniref90_annotation.tsv")
	mymsp_taxa = os.path.join(main_folder, config.basename + "_MSPminer_msp.taxonomy.tsv")
	mymsp_ann = os.path.join(main_folder, config.basename + "_MSPminer_annotation.tsv")
	mymsp_ann_taxa = os.path.join(main_folder, config.basename + "_MSPminer_annotation.taxonomy.tsv")
	mymsp_detail_family = os.path.join(main_folder, config.basename + "_MSPminer_proteinfamilies.detail.tsv")
	mymsp_detail = os.path.join(main_folder, config.basename + "_MSPminer_proteinfamilies.ORF.detail.tsv")
	msp_detail_family = os.path.join(main_folder, config.basename + "_MSP_proteinfamilies.detail.tsv")
	msp_detail = os.path.join(main_folder, config.basename + "_MSP_proteinfamilies.ORF.detail.tsv")
	mymsp_detail_taxa_family = os.path.join(main_folder, config.basename + "_proteinfamilies_annotation.taxonomy.tsv")
	mymsp_detail_taxa = os.path.join(main_folder, config.basename + "_protein_annotation.taxonomy.tsv")

	gene_rpk = os.path.join(main_folder, config.basename + "_genecatalogs_counts.RPK.tsv")
	gene_relab = os.path.join(main_folder, config.basename + "_genecatalogs_nrm.tsv")
	family_count_all = os.path.join(main_folder, config.basename + "_proteinfamilies_counts.all.tsv")
	family_count = os.path.join(main_folder, config.basename + "_proteinfamilies_counts.tsv")
	family_rpk = os.path.join(main_folder, config.basename + "_proteinfamilies_counts.RPK.tsv")
	family_relab = os.path.join(main_folder, config.basename + "_proteinfamilies_nrm.tsv")
	abundance_ann_family =  os.path.join(main_folder, config.basename + "_DNA_proteinfamilies.abundance.detail.tsv")
	abundance_ann = os.path.join(main_folder, config.basename + "_DNA_proteinfamilies.ORF.abundance.detail.tsv")

	DA = os.path.join(main_folder, "DA")
	family_smooth = os.path.join(DA, config.basename + "_proteinfamilies_nrm.smooth.tsv")
	DA_metadata = os.path.join(DA, "maaslin2_metadata.tsv")
	DA_results = os.path.join(config.maaslin2_dir, "all_results.tsv")
	DA_prefix = os.path.join(DA, config.basename + "_stat_diff_family_abundance.tsv")
	stat_results = os.path.join(DA, config.basename + "_stat_diff_family_abundance.stat.tsv")
	stat_results_minq = os.path.join(DA, config.basename + "_stat_diff_family_abundance.minQ.stat.tsv")
	fold_results = os.path.join(DA, config.basename + "_stat_diff_family_abundance.fold.tsv")
	fold_results_minq = os.path.join(DA, config.basename + "_stat_diff_family_abundance.minQ.fold.tsv")
	pre_results = os.path.join(DA, config.basename + "_stat_diff_family_abundance.prevalence.tsv")
	pre_results_minq = os.path.join(DA, config.basename + "_stat_diff_family_abundance.minQ.prevalence.tsv")
	summary_stat = os.path.join(DA, config.basename + "_diff_family_abundance.stat.tsv")
	summary_stat_minq = os.path.join(DA, config.basename + "_diff_family_abundance.minQ.stat.tsv")
	DA_detail =  os.path.join(main_folder, config.basename + "_MaAsLin2_proteinfamilies.DA.detail.tsv")
	DA_detail_minq = os.path.join(main_folder, config.basename + "_MaAsLin2_proteinfamilies.DA-minQ.detail.tsv")


	myname = os.path.basename(taxonomy_annotation)
	taxa_family = os.path.join(main_folder, myname)
	taxa = re.sub("proteinfamilies", "protein", taxa_family)
	myprotein_family_ann = {}
	myprotein_ann = {}
	time_equation = config.time  # xxx hours defined in global config
	mem_equation = config.memory  # xxx GB defined in global config


	#### MSP annotation ####
	if abundance_conf["mspminer"] != "no" and abundance_conf["mspminer"] != "No":
		## summary gene catalog abundance ##
		# remove eukaryotic contamination
		mylog = re.sub(".tsv", ".log", count_file)
		uniref_taxonomy_tmp = os.path.join(tmps_dir, os.path.basename(uniref_taxonomy))
		workflow.add_task_gridable(
				"metawibele_abundance_filtering -a [depends[0]] -f good -i [depends[1]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [gene_catalog_count, uniref_taxonomy_tmp, TrackedExecutable("metawibele_abundance_filtering")],
				targets = [count_file],
				args = [mylog],
				cores = 1,
				time = time_equation,
				mem = mem_equation,
				name = "abundance_filtering")
		
		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
	            depends = [count_file],
	            targets = [os.path.join(tmps_dir, os.path.basename(count_file))],
		        cores = 1,
	            name = "ln__counts_file")
		
		## binning co-abundant genes across metagenomic samples ##
		# run MSPminer
		raw_conf = abundance_conf["mspminer"]
		myconf = os.path.join(main_folder, os.path.basename(raw_conf))
		#myconf = re.sub(".cfg", ".refined.cfg", raw_conf)
		msp_output = main_folder
		open_file = open(raw_conf, "r")
		open_out = open(myconf, "w")
		for line in open_file.readlines():
			line = line.strip()
			if not len(line):
				open_out.write(line + "\n")
				continue
			if re.search("output_dir", line):
				mym = re.search("output_dir\=([\S]+)", line)
				myraw = mym.group(1)
				msp_output = os.path.join(msp_output, myraw)
				line = re.sub(myraw, msp_output, line)
				os.system("rm -rf " + msp_output)
				os.system("mkdir -p " + msp_output)
			if re.search("count_matrix_file\=", line):
				line = "count_matrix_file=" + count_file
			open_out.write(line + "\n")
		# foreach line
		open_file.close()
		open_out.close()

		mylog = re.sub(".tsv", ".run_mspminer.log", mymsp)
		myoutfile = os.path.join(msp_output, "all_seeds.tsv")
		workflow.add_task(
				"mspminer [depends[0]] >[args[0]] 2>&1",
				depends = [myconf, count_file, TrackedExecutable("mspminer")],
				targets = [myoutfile],
				args = [mylog],
				cores = threads,
				name = "mspminer")

		# summary MSPs
		mylog = re.sub(".tsv", ".log", mymsp)
		workflow.add_task(
				"metawibele_mspminer_msp -a [depends[0]] -p [args[0]] -o [targets[0]] >[args[1]] 2>&1",
				depends = [uniref_taxonomy, myoutfile, TrackedExecutable("metawibele_mspminer_msp")],
				targets = [mymsp],
				args = [msp_output, mylog],
				cores = 1,
				name = "mspminer_msp")

		## taxonomy annotation for MSPs ##
		# summary taxon for each gene
		mylog = re.sub(".tsv", ".log", mymsp_uniref)
		workflow.add_task(
				"metawibele_mspminer_msp_uniref_annotation -a [depends[0]] -i [depends[1]] -o [targets[0]] >[args[0]] 2>&1",
				depends = [uniref_taxonomy, mymsp, TrackedExecutable("metawibele_mspminer_msp_uniref_annotation")],
				targets = [mymsp_uniref],
				args = [mylog],
				cores = 1,
				name = "mspminer_msp_uniref_annotation")

		# annotate taxon for each MSP
		mylog = re.sub(".tsv", ".log", mymsp_taxa)
		workflow.add_task(
			"metawibele_mspminer_msp_taxonomy_annotation -a [depends[0]] -t UniRef90_homology -o [targets[0]] >[args[0]] 2>&1",
			depends = [mymsp_uniref, TrackedExecutable("metawibele_mspminer_msp_taxonomy_annotation")],
			targets = [mymsp_taxa],
			args = [mylog],
			cores = 1,
			name = "mspminer_msp_taxonomy_annotation")

		# summary MSP annotation for each gene
		mylog = re.sub(".tsv", ".log", mymsp_ann)
		workflow.add_task(
			"metawibele_mspminer_protein -a [depends[0]] -m [depends[1]] -g [depends[2]] -o [targets[0]] >[args[0]] 2>&1",
			depends = [mymsp_taxa, mymsp, uniref_taxonomy, TrackedExecutable("metawibele_mspminer_protein")],
			targets = [mymsp_ann, mymsp_ann_taxa],
			args = [mylog],
			cores = 1,
			name = "mspminer_protein")

		# summary annotation for protein family
		mylog = re.sub(".tsv", ".log", mymsp_detail_family)
		workflow.add_task(
			"metawibele_mspminer_protein_family -f centroid -m [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
			depends = [mymsp_ann, protein_family, TrackedExecutable("metawibele_mspminer_protein_family")],
			targets = [mymsp_detail_family, mymsp_detail],
			args = [mylog],
			cores = 1,
			name = "mspminer_protein_family")
		
		mylog = re.sub(".tsv", ".log", msp_detail_family)
		workflow.add_task(
			"metawibele_msp_protein_family -m [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
			depends = [mymsp_ann_taxa, protein_family, TrackedExecutable("metawibele_msp_protein_family")],
			targets = [msp_detail_family, msp_detail],
			args = [mylog],
			cores = 1,
			name = "msp_protein_family")

		mylog = re.sub(".tsv", ".log", mymsp_detail_taxa_family)
		workflow.add_task(
			"metawibele_mspminer_protein_family_taxonomy -a [depends[0]] -s [args[0]] -o [targets[0]] >[args[1]] 2>&1",
			depends = [mymsp_ann_taxa, protein_family, TrackedExecutable("metawibele_mspminer_protein_family_taxonomy")],
			targets = [mymsp_detail_taxa_family, mymsp_detail_taxa],
			args = [config.taxa_final, mylog],
			cores = 1,
			name = "mspminer_protein_family_taxonomy")

		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
		    depends=[mymsp_detail_taxa_family],
			targets=[taxonomy_annotation_family],
		    cores = 1,
		    name="ln__mymsp_detail_taxa_family")
		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
			depends=[mymsp_detail_taxa],
		    targets=[taxonomy_annotation],
		    cores=1,
		    name="ln__mymsp_detail_taxa")

		myprotein_family_ann[mymsp_detail_family] = ""
		myprotein_ann[mymsp_detail] = ""
	# if MSP

	## DNA abundance ##
	if abundance_conf["dna_abundance"] != "no" and abundance_conf["dna_abundance"] != "No":

		# get raw counts for protein families
		mylog = re.sub(".tsv", ".log", family_count_all)
		workflow.add_task_gridable(
			"metawibele_sum_to_protein_family_abundance -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1",
			depends = [count_file, protein_family, TrackedExecutable("metawibele_sum_to_protein_family_abundance")],
			targets = [family_count_all],
			args = [mylog],
			cores = 1,
			time = time_equation,
			mem = mem_equation,
			name = "sum_to_protein_family_abundance")

		# filter out bad protein
		mylog = re.sub(".tsv", ".log", family_count)
		uniref_taxonomy_family_tmp = os.path.join(tmps_dir, os.path.basename(uniref_taxonomy_family))
		workflow.add_task(
			"metawibele_abundance_filtering -i [depends[0]] -f good -a [depends[1]] -o [targets[0]] >[args[0]] 2>&1",
			depends = [uniref_taxonomy_family_tmp, family_count_all, TrackedExecutable("metawibele_abundance_filtering")],
			targets = [family_count],
			args = [mylog],
			cores = 1,
			name = "abundance_filtering")

		# normalized to RPK and relative abundance
		mylog = re.sub(".tsv", ".log", family_rpk)
		workflow.add_task(
			"metawibele_abundance_RPK -i [depends[0]] -l [depends[1]] -t cluster -o [targets[0]] >[args[0]] 2>&1",
			depends = [family_count, protein_family, TrackedExecutable("metawibele_abundance_RPK")],
			targets = [family_rpk],
			args = [mylog],
			cores = 1,
			name = "abundance_RPK")

		mylog = re.sub(".tsv", ".log", family_relab)
		workflow.add_task(
			"metawibele_abundance_normalization -i [depends[0]] -u [args[0]] -o [targets[0]] >[args[1]] 2>&1",
			depends = [family_rpk, TrackedExecutable("metawibele_abundance_normalization")],
			targets = [family_relab],
			args = [config.normalize, mylog],
			cores = 1,
			name = "abundance_normalization")

		mylog = re.sub(".tsv", ".log", gene_rpk)
		workflow.add_task_gridable(
			"metawibele_abundance_RPK -i [depends[0]] -l [depends[1]] -t fasta -o [targets[0]] >[args[0]] 2>&1",
			depends = [count_file, gene_catalog_seq, TrackedExecutable("metawibele_abundance_RPK")],
			targets = [gene_rpk],
			args = [mylog],
			cores = 1,
			time = time_equation,
			mem = mem_equation,
			name = "abundance_RPK_gene")

		mylog = re.sub(".tsv", ".log", gene_relab)
		workflow.add_task_gridable(
			"metawibele_abundance_normalization -i [depends[0]] -u [args[0]] -o [targets[0]] >[args[1]] 2>&1",
			depends = [gene_rpk, TrackedExecutable("metawibele_abundance_normalization")],
			targets = [gene_relab],
			args = [config.normalize, mylog],
			cores = 1,
			time = time_equation,
			mem = mem_equation,
			name = "abundance_normalization")

		# format abundance annotation
		mylog = re.sub(".tsv", ".log", abundance_ann_family)
		workflow.add_task(
			"metawibele_abundance_annotator -a [depends[0]] -f protein_family -t [args[0]] -o [targets[0]] >[args[1]] 2>&1",
			depends = [family_relab, TrackedExecutable("metawibele_abundance_annotator")],
			targets = [abundance_ann_family],
			args = [abundance_conf["dna_abundance"], mylog],
			cores = 1,
			name = "abundance_annotator")

		mylog = re.sub(".tsv", ".log", abundance_ann)
		workflow.add_task_gridable(
			"metawibele_abundance_annotator -a [depends[0]] -f protein -t [args[0]] -o [targets[0]] >[args[1]] 2>&1",
			depends = [gene_relab, TrackedExecutable("metawibele_abundance_annotator")],
			targets = [abundance_ann],
			args = [abundance_conf["dna_abundance"], mylog],
			cores = 1,
			time = time_equation,
			mem = mem_equation,
			name = "abundance_annotator")

		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
			depends = [family_relab],
		    targets = [protein_family_relab],
		    cores = 1,
		    name="ln__proteinfamilies_relab")
		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
			depends = [family_relab],
		    targets = [os.path.join(tmps_dir, os.path.basename(family_relab))],
		    cores = 1,
		    name="ln__proteinfamilies_relab")

		myprotein_family_ann[abundance_ann_family] = ""
		myprotein_ann[abundance_ann] = ""
	# if DNA abundance


	## Differential abundance for DNA abundance##
	if abundance_conf["dna_da"] != "no" and abundance_conf["dna_da"] != "No":
		os.system("mkdir -p " + DA)

		## prepare data for maaslin2 ##
		mylog = re.sub(".tsv", ".log", family_smooth)
		family_relab_tmp = os.path.join(tmps_dir, os.path.basename(family_relab))
		workflow.add_task_gridable(
				"metawibele_abundance_smoothing -i [depends[0]] -f 0 -o [targets[0]] >[args[1]] 2>&1",
				depends = [family_relab_tmp, TrackedExecutable("metawibele_abundance_smoothing")],
				targets = [family_smooth],
				args = [config.tshld_prevalence, mylog],
				cores = 1,
				time = time_equation,
				mem = mem_equation,
				name = "abundance_smoothing")
		
		feature_pcl = re.sub(".tsv", ".feature.pcl", family_smooth)
		workflow.add_task("ln -fs [depends[0]] [targets[0]]",
		        depends = [family_smooth],
		        targets = [feature_pcl],
		        cores = 1,
		        name="ln__proteinfamilies_relab_feature")

		feature_tsv = re.sub(".pcl", ".tsv", feature_pcl)
		mylog = re.sub(".pcl", ".pcl.log", feature_tsv)
		workflow.add_task(
				"metawibele_transpose < [depends[0]] > [targets[0]] 2> [args[0]]",
				depends = [feature_pcl, TrackedExecutable("metawibele_transpose")],
				targets = [feature_tsv],
				args = [mylog],
				cores = 1,
				name = "transpose")

		mylog = re.sub(".tsv", ".log", DA_metadata)
		workflow.add_task(
				"metawibele_metadata_format -i [depends[0]] -o [targets[0]] >[args[0]] 2>&1 ",
				depends = [config.metadata, TrackedExecutable("metawibele_metadata_format")],
				targets = [DA_metadata],
				args = [mylog],
				cores = 1,
				name = "metadata_format")

		## run DA  ##
		mylog = re.sub(".tsv", ".results.log", DA_metadata)
		myresults = re.sub(".tsv", ".fdr_correction.correct_per_level.tsv", DA_results)
		workflow.add_task(
				"metawibele_maaslin2 -i [depends[0]] -m [depends[1]] -n [args[2]] -w [args[0]] -o [args[1]] > [args[3]] 2>&1",
				depends = [feature_tsv, DA_metadata, TrackedExecutable("metawibele_maaslin2")],
				targets = [myresults],
				args = [os.path.join(DA, "maaslin2_output"), "all_results.tsv", threads, mylog],
				cores = threads,
				name = "maaslin2")

		## collect results info ##
		mylog = re.sub(".tsv", ".tsv.log", stat_results)
		workflow.add_task(
				"metawibele_maaslin2_collection -a [depends[0]] -s [depends[1]] -i [depends[2]] -o [args[0]] > [args[1]] 2>&1",
				depends = [family_relab, family_smooth, myresults, TrackedExecutable("metawibele_maaslin2_collection")],
				targets = [stat_results, fold_results, pre_results, stat_results_minq, fold_results_minq, pre_results_minq],
				args = [DA_prefix, mylog],
				cores = 1,
				name = "maaslin2_collection")

		## summary DA results ##
		myeffect = config.effect_size
		myeffect = re.sub("\(", "_", myeffect)
		myeffect = re.sub("\)", "", myeffect)
		mylog = re.sub(".tsv", ".tsv.log", summary_stat)
		summary_stat_all = re.sub(".tsv", ".all.tsv", summary_stat)
		workflow.add_task(
				"metawibele_maaslin2_summary -a [depends[0]] -b [depends[1]] -c [depends[2]] -p [args[0]] -q [args[1]] -e [args[2]] -o [targets[0]] > [args[3]] 2>&1",
				depends = [stat_results, fold_results, pre_results, TrackedExecutable("metawibele_maaslin2_summary")],
				targets = [summary_stat, summary_stat_all],
				args = [config.tshld_prevalence, config.tshld_qvalue, myeffect, mylog],
				cores = 1,
				name = "maaslin2_summary")
		
		mylog = re.sub(".tsv", ".tsv.log", summary_stat_minq)
		summary_stat_minq_all = re.sub(".tsv", ".all.tsv", summary_stat_minq)
		workflow.add_task(
				"metawibele_maaslin2_summary -a [depends[0]] -b [depends[1]] -c [depends[2]] -p [args[0]] -q [args[1]] -e [args[2]] -o [targets[0]] > [args[3]] 2>&1",
				depends = [stat_results_minq, fold_results_minq, pre_results_minq, TrackedExecutable("metawibele_maaslin2_summary")],
				targets = [summary_stat_minq, summary_stat_minq_all],
				args = [config.tshld_prevalence, config.tshld_qvalue, myeffect, mylog],
				cores = 1,
				name = "maaslin2_summary")

		mytag = abundance_conf["dna_da"] + "-sig"
		DA_detail_sig = re.sub(".detail.tsv", "-sig.detail.tsv", DA_detail)
		mylog = re.sub(".tsv", ".tsv.log", DA_detail_sig)
		workflow.add_task(
				"metawibele_maaslin2_annotator -s [depends[0]] -t [args[0]] -e [args[1]] -o [targets[0]] > [args[2]] 2>&1",
				depends = [summary_stat, TrackedExecutable("metawibele_maaslin2_annotator")],
				targets = [DA_detail_sig],
				args = [mytag, myeffect, mylog],
				cores = 1,
				name = "maaslin2_annotator")
		myprotein_family_ann[DA_detail_sig] = ""
		#summary_stat_all = re.sub(".tsv", ".all.tsv", summary_stat)
		mylog = re.sub(".tsv", ".tsv.log", DA_detail)
		workflow.add_task(
				"metawibele_maaslin2_annotator -s [depends[0]] -t [args[0]] -e [args[1]] -o [targets[0]] > [args[2]] 2>&1",
				depends = [summary_stat_all, TrackedExecutable("metawibele_maaslin2_annotator")],
				targets = [DA_detail],
				args = [abundance_conf["dna_da"], myeffect, mylog],
				cores = 1,
				name = "maaslin2_annotator")
		myprotein_family_ann[DA_detail] = ""

		mytag = abundance_conf["dna_da"] + "-minQ-sig"
		DA_detail_minq_sig = re.sub(".detail.tsv", "-sig.detail.tsv", DA_detail_minq)
		mylog = re.sub(".tsv", ".tsv.log", DA_detail_minq_sig)
		workflow.add_task(
				"metawibele_maaslin2_annotator -s [depends[0]] -t [args[0]] -e [args[1]] -o [targets[0]] > [args[2]] 2>&1",
				depends = [summary_stat_minq, TrackedExecutable("metawibele_maaslin2_annotator")],
				targets = [DA_detail_minq_sig],
				args = [mytag, myeffect, mylog],
				cores = 1,
				name = "maaslin2_annotator")
		myprotein_family_ann[DA_detail_minq_sig] = ""
		mytag = abundance_conf["dna_da"] + "-minQ"
		#summary_stat_minq_all = re.sub(".tsv", ".all.tsv", summary_stat_minq)
		mylog = re.sub(".tsv", ".tsv.log", DA_detail_minq)
		workflow.add_task(
				"metawibele_maaslin2_annotator -s [depends[0]] -t [args[0]] -e [args[1]] -o [targets[0]] > [args[2]] 2>&1",
				depends = [summary_stat_minq_all, TrackedExecutable("metawibele_maaslin2_annotator")],
				targets = [DA_detail_minq],
				args = [mytag, myeffect, mylog],
				cores = 1,
				name = "maaslin2_annotator")
		myprotein_family_ann[DA_detail_minq] = ""

	# if DA for DNA abundance


	if len(myprotein_family_ann.keys()) > 0:
		for myfile in myprotein_family_ann.keys():
			protein_family_ann_list[myfile] = ""
	if len(myprotein_ann.keys()) > 0:
		for myfile in myprotein_ann.keys():
			protein_ann_list[myfile] = ""

	return myprotein_family_ann, myprotein_ann, main_folder


def integration_annotation (workflow, integration_conf,
                          protein_family_ann_list, protein_ann_list,
                          uniref_annotation_family, uniref_annotation,
                          taxonomy_annotation_family, taxonomy_annotation,
                          threads, output_folder, protein_family_ann, protein_family_attr):

	"""
	Finalization annotation workflow

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		integration_conf: configuration for integration annotations
		protein_family_ann_list_file: List file of annotations files for protein families
		protein_ann_list_file: List file of annotations files for protein
		uniref_annotation_family: uniref annotation for protein families
		uniref_annotation: uniref annotation for protein
		taxonomy_annotation_family: taxonomic annotation for protein families
		taxonomy_annotation: taxonomic annotation for protein families
		threads (int): The number of threads/cores for clustering to use.
		output_folder (string): The path of the output folder.
		protein_family_ann: integration annotation file
		protein_family_attr: integration annotation attribution file

	Requires:
		protein_family_ann_list_file: List file of annotations files for protein families
		protein_ann_list_file: List file of annotations files for protein
		uniref_annotation_family: uniref annotation for protein families
		uniref_annotation: uniref annotation for protein
		taxonomy_annotation_family: taxonomic annotation for protein families
		taxonomy_annotation: taxonomic annotation for protein families

	Returns:
		files: integration annotation file

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		protein_family_ann, protein_family_attr, annotation_output_folder = characterization.integration_annotation (workflow,
	                                                                                                          protein_family_ann_list_file, protein_ann_list_file,
	                                                                                                          uniref_annotation_family, uniref_annotation,
	                                                                                                          taxonomy_annotation_family, taxonomy_annotation,
	                                                                                                          args.threads,
	                                                                                                          output_dir, protein_family_ann, protein_family_attr)

		# run the workflow
		workflow.go()

	"""

	# define the annotation output files
	main_folder = output_folder
	summary_ann_family = os.path.join(main_folder, config.basename + "_summary_proteinfamilies_annotation.tsv")
	summary_ann = os.path.join(main_folder, config.basename + "_summary_protein_annotation.tsv")
	all_ann_family = os.path.join(main_folder, config.basename + "_summary_proteinfamilies_annotation.all.tsv")
	all_ann = os.path.join(main_folder, config.basename + "_summary_protein_annotation.all.tsv")
	final_ann_family = os.path.join(main_folder, config.basename + "_proteinfamilies_annotation.tsv")
	final_attr_family = os.path.join(main_folder, config.basename + "_proteinfamilies_annotation.attribute.tsv")
	final_ann = os.path.join(main_folder, config.basename + "_protein_annotation.tsv")
	final_attr = os.path.join(main_folder, config.basename + "_protein_annotation.attribute.tsv")

	# collect list
	proteinfamily_list = []
	protein_list = []
	protein_family_ann_list_file = os.path.join(output_folder, "protein_family_ann.list")
	protein_ann_list_file = os.path.join(output_folder, "protein_ann.list")
	utilities.dict_to_file (protein_family_ann_list, protein_family_ann_list_file)
	utilities.dict_to_file (protein_ann_list, protein_ann_list_file)
	proteinfamily_list = protein_family_ann_list
	protein_list = protein_ann_list
	time_equation = config.time  # xxx hours defined in global config
	mem_equation = config.memory  # xxx GB defined in global config

	# summarize annotationa
	if integration_conf["summary_ann"].lower() == "yes":
		myinputs = []
		myinputs.append(protein_family_ann_list_file)
		myinputs.append(uniref_annotation_family)
		myinputs.extend(proteinfamily_list)
		mylog = re.sub(".tsv", ".tsv.log", summary_ann_family)
		workflow.add_task(
				"metawibele_summary_function_annotation -l [depends[0]] -a [depends[1]] -o [targets[0]] > [args[0]] 2>&1",
				depends = utilities.add_to_list(myinputs, TrackedExecutable("metawibele_summary_function_annotation")),
				targets = [summary_ann_family],
				args = [mylog],
				cores = 1,
				name = "summary_function_annotation")

		mylog = re.sub(".tsv", ".tsv.log", summary_ann)
		myinputs = []
		myinputs.append(protein_ann_list_file)
		myinputs.append(uniref_annotation)
		myinputs.extend(protein_list)
		workflow.add_task(
				"metawibele_summary_function_annotation -l [depends[0]] -a [depends[1]] -o [targets[0]] > [args[0]] 2>&1",
				depends = utilities.add_to_list(myinputs, TrackedExecutable("metawibele_summary_function_annotation")),
				targets = [summary_ann],
				args = [mylog],
				cores = 1,
				name = "summary_function_annotation")

		mylog = re.sub(".tsv", ".tsv.log", all_ann_family)
		workflow.add_task(
				"metawibele_summary_all_annotation -a [depends[0]] -t [depends[1]] -o [targets[0]] > [args[0]] 2>&1",
				depends = [summary_ann_family, taxonomy_annotation_family, TrackedExecutable("metawibele_summary_all_annotation")],
				targets = [all_ann_family],
				args = [mylog],
				cores = 1,
				name = "summary_all_annotation")

		mylog = re.sub(".tsv", ".tsv.log", all_ann)
		workflow.add_task(
				"metawibele_summary_all_annotation -a [depends[0]] -t [depends[1]] -o [targets[0]] > [args[0]] 2>&1",
				depends = [summary_ann, taxonomy_annotation, TrackedExecutable("metawibele_summary_all_annotation")],
				targets = [all_ann],
				args = [mylog],
				cores = 1,
				name = "summary_all_annotation")

	## finalize annotation
	if integration_conf["finalization"].lower() == "yes": 
		mylog = re.sub(".tsv", ".tsv.log", final_ann_family)
		myinputs = []
		myinputs.append(protein_family_ann_list_file)
		myinputs.append(summary_ann_family)
		myinputs.append(taxonomy_annotation_family)
		myinputs.append(uniref_annotation_family)
		myinputs.extend(proteinfamily_list)
		workflow.add_task(
				"metawibele_finalize_annotation -l [depends[0]] -a [depends[1]] -t [depends[2]] -u [depends[3]] -s protein_family -o [targets[0]] > [args[0]] 2>&1",
				depends = utilities.add_to_list(myinputs, TrackedExecutable("metawibele_finalize_annotation")),
				targets = [final_ann_family, final_attr_family],
				args = [mylog],
				cores = 1,
				name = "finalize_annotation")

		mylog = re.sub(".tsv", ".tsv.log", final_ann)
		myinputs = []
		myinputs.append(protein_ann_list_file)
		myinputs.append(summary_ann)
		myinputs.append(taxonomy_annotation)
		myinputs.append(uniref_annotation)
		myinputs.extend(protein_list)
		workflow.add_task(
				"metawibele_finalize_annotation -l [depends[0]] -a [depends[1]] -t [depends[2]] -u [depends[3]] -s protein -o [targets[0]] > [args[0]] 2>&1",
				depends = utilities.add_to_list(myinputs, TrackedExecutable("metawibele_finalize_annotation")),
				targets = [final_ann, final_attr],
				args = [mylog],
				cores = 1,
				name = "finalize_annotation")

	return final_ann_family, final_attr_family, main_folder

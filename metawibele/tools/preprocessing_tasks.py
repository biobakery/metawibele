#!/usr/bin/env python3

"""
Preprocessing Workflow: assembly metagenomic shotgun sequencing reads and build gene catalogs
A collection of tasks for workflows with gene catalogs

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
from metawibele import utilities, config, files


def assembly (workflow, input_dir, sample_file, extension_paired, extension_orphan,
             threads, output_folder, contigs):
	"""
	This set of tasks will run assembly on the input files provided.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		input_dir: The direcory path of fastq files.
		sample_file: The sample list file.
		extension_paired: The extension for paired reads.
		extension_orphan: The extension for arphan reads.
		threads (int): The number of threads/cores for clustering to use.
		output_folder (string): The path of the output folder.
		contigs: The summarized contig file.

	Requires:
		metahit v1.1.3: A program for assembling metagenomic sequencing reads
		fastq files

	Returns:
		string: the name of contigs file.

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add assembly tasks
		mycontigs  = preprocessing_tasks.assembly (workflow, input_dir, args.sample_file,
												   args.extension_paired, args.extension_orphan,
												   args.threads,
												   assembly_dir, contigs)
		# run the workflow
		workflow.go()
	"""

	# ================================================
	# collect sequences
	# ================================================
	samples = []
	open_file = open(sample_file, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		samples.append(info[0])
	# foreach sample
	open_file.close()
	split_dir = input_dir
	assembly_dir = output_folder

	split_files = []
	contigs_list = []
	tmp1 = extension_paired.split(",")
	tmp2 = extension_orphan.split(",")
	for sample in samples:
		#f_seq = os.path.join(split_dir, sample + tmp1[0])
		#r_seq = os.path.join(split_dir, sample + tmp1[1])
		mypair = "none"
		mypair_tmp = []
		for item in tmp1:
			if item == "none":
				continue
			mypair_tmp.append(os.path.join(split_dir, sample + item))
		if len(mypair_tmp) > 0:
			mypair_tmp = utilities.paired_reads(mypair_tmp, "fastq.gz", pair_identifier=".R1")
			mypair = ",".join(mypair_tmp)	
		myorphan = "none"
		for item in tmp2:
			if item == "none":
				continue
			if myorphan == "none":
				myorphan = os.path.join(split_dir, sample + item)
			else:
				myorphan = myorphan + "," + os.path.join(split_dir, sample + item)
		#split_files.append((sample, f_seq, r_seq, myorphan))
		split_files.append((sample, mypair, myorphan))
		
		seq_base = sample
		megahit_contig_dir = os.path.join(assembly_dir, seq_base)
		megahit_contig = os.path.join(megahit_contig_dir, '%s.contigs.fa' % seq_base)
		contigs_list.append(megahit_contig)

	## run MEGAHIT
	os.system("mkdir -p " + assembly_dir)
	for (sample, mypair, myorphan) in split_files:
		seq_base = sample
		megahit_contig_dir = os.path.join(assembly_dir, seq_base)
		megahit_contig = os.path.join(megahit_contig_dir, '%s.contigs.fa' % seq_base)

		## MEGAHIT needs memory in a byte format so let's take care of data
		time_equation = "24*60 if file_size('[depends[0]]') < 25 else 6*24*60" # 24 hours or more depending on file size
		mem_equation = "32*1024 if file_size('[depends[0]]') < 25 else 3*32*1024" # 32 GB or more depending on file size
		mylog = os.path.join(assembly_dir, '%s.log' % seq_base)
		tmp = mypair.split(",")
	
		if len(tmp) == 2:	# paired reads:		
			f_seq = tmp[0]
			r_seq = tmp[1]
			if myorphan != "none":
				workflow.add_task_gridable("megahit -1 [depends[0]] -2 [depends[1]] -r [args[2]] -t [args[0]] -o [args[3]] --out-prefix [args[1]] >[args[4]] 2>&1",
									depends = [f_seq, r_seq, TrackedExecutable("megahit")],
									targets = [megahit_contig],
									args = [threads, seq_base, myorphan, megahit_contig_dir, mylog],
									cores = threads,
									mem = mem_equation,
									time = time_equation,
									name = sample + "__megahit")
			else:
				workflow.add_task_gridable("megahit -1 [depends[0]] -2 [depends[1]] -t [args[0]] -o [args[2]] --out-prefix [args[1]] >[args[3]] 2>&1",
									depends = [f_seq, r_seq, TrackedExecutable("megahit")],
									targets = [megahit_contig],
									args = [threads, seq_base, megahit_contig_dir, mylog],
									cores = threads,
									name = sample + "__megahit")
		else:
			workflow.add_task_gridable("megahit -r [depends[0]] -t [args[0]] -o [args[2]] --out-prefix [args[1]] >[args[3]] 2>&1",
									depends = [mypair, TrackedExecutable("megahit")],
									targets = [megahit_contig],
									args = [threads, seq_base, megahit_contig_dir, mylog],
									cores = threads,
									name = sample + "__megahit")
			
		
	for myfile in contigs_list:
		mym = re.search("([^\/]+)$", myfile)
		myname = mym.group(1)
		myfile_new = assembly_dir + "/" + myname
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [myfile],
			targets = [myfile_new],
			name = "ln__" + myname)


	## combine contigs sequences
	mylog = re.sub(".tsv", ".log", contigs)
	workflow.add_task(
		"format_contig_sequences.py -p [args[0]] -e contigs.fa -o [targets[0]] > [args[1]] 2>&1",
		depends=contigs_list,
		targets=[contigs],
		args=[assembly_dir, mylog],
		name="format_contig_table")

	return contigs_list


def gene_calling (workflow, assembly_dir, assembly_extentsion, sample_file,
                 prokka_dir, prodigal_dir,
                 threads,
                 gene_file, gene_PC_file, protein_file, protein_sort,
                 gene_info, complete_gene, complete_protein):
	"""
	This set of tasks will run gene-calling workflow.

	Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		assembly_dir: The direcory path of assembly results.
		sample_file: The sample list file.
		prokka_dir: The direcory path of prokka results.
		prodigal_dir: The direcory path of prodigal results.
		gene_file: The fasta file of gene nucleotide sequences.
		gene_PC_file: The fasta file of protein coding gene nucleotide sequences.
		protein_file: The fasta file of protein sequences.
		protein_sort: The sorted fasta file of protein sequences.
		gene_info: The summaized gene calling file.
		complete_gene: The fasta file of gene nucleotide sequences for complete ORFs.
		complete_protein: The fasta file of protein sequences for complete ORFs.

	Requires:
		prokka 1.14-dev: rapid prokaryotic genome annotation
		prodigal v2.6: gene prediction
		usearch (tested with usearch v9.0.2132_i86linux64)
		assembled contig files

	Returns:
		string: name of gene files

	Example:
		from anadama2 import Workflow
		from MetaWIBELE.characterize import characterization

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add gene calling tasks
		mygene, myprotein = preprocessing_tasks.gene_calling (workflow, assembly_dir, args.sample_file,
															  prokka_dir, prodigal_dir,
															  gene_file, gene_PC_file, protein_file, protein_sort,
															  gene_info, complete_gene, complete_protein)
		# run the workflow
		workflow.go()
	"""

	# ================================================
	# collect sequences
	# ================================================
	samples = []
	sequence_files = []
	open_file = open(sample_file, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myfile = os.path.join(assembly_dir, info[0], info[0] + "%s" % assembly_extentsion)
		samples.append(info[0])
		sequence_files.append(myfile)
	# foreach sample
	open_file.close()

	filtered_contigs = sequence_files

	# ================================================
	# Gene calling
	# ================================================
	os.system("mkdir -p " + prodigal_dir)
	time_equation = "24*60 if file_size('[depends[0]]') < 25 else 6*24*60" # 24 hours or more depending on file size
	mem_equation = "32*1024 if file_size('[depends[0]]') < 25 else 3*32*1024" # 32 GB or more depending on file size
	fna_file = []
	faa_file = []
	fna_file_tmp = []
	faa_file_tmp = []

	## Using Prodigal
	for contig in filtered_contigs:
		contig_base = os.path.basename(contig).split(os.extsep)[0]
		annotation_dir = prodigal_dir + "/" + contig_base
		os.system("mkdir -p " + annotation_dir)
		gff_file = os.path.join(annotation_dir, '%s.gff' % contig_base)
		cds_file = os.path.join(annotation_dir, '%s.fna' % contig_base)
		cds_aa = os.path.join(annotation_dir, '%s.faa' % contig_base)
		score = os.path.join(annotation_dir, '%s.gene_score.txt' % contig_base)
		stdout_log = os.path.join(annotation_dir, '%s.stdout.log' % contig_base)
		faa_file_tmp.append(cds_aa)

		workflow.add_task_gridable('prodigal -m -p meta -i [depends[0]] '
		                           '-f gff -o [targets[0]] -d [targets[1]] -s [targets[3]] '
		                           '-a [targets[2]] '
		                           '>[args[0]] 2>&1',
								   depends = [contig, TrackedExecutable("prodigal")],
								   targets = [gff_file, cds_file, cds_aa, score],
								   args = [stdout_log],
								   cores = threads,
								   mem = mem_equation,
								   time = time_equation,
								   name = contig_base + "__prodigal")
	
	for myfile in faa_file_tmp:
		mym = re.search("([^\/]+)$", myfile)
		myname = mym.group(1)
		myfile_new = prodigal_dir + "/" + myname
		faa_file.append(myfile_new)
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [myfile],
			targets = [myfile_new],
			name = "ln__" + myname)
		myfna = re.sub(".faa", ".fna", myfile)
		myfna_new = re.sub(".faa", ".fna", myfile_new)
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [myfna],
			targets = [myfna_new],
			name = "ln__" + myname)
		mygff = re.sub(".faa", ".gff", myfile)
		mygff_new = re.sub(".faa", ".gff", myfile_new)
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [mygff],
			targets = [mygff_new],
			name = "ln__" + myname)



	## Calling genes with Prokka
	os.system("mkdir -p " + prokka_dir)

	for contig in filtered_contigs:
		contig_base = os.path.basename(contig).split(os.extsep)[0]
		mym = re.search("([^\/]+$)", contig_base)
		sample = mym.group(1)
		annotation_dir = prokka_dir + "/" + sample
		os.system("mkdir -p " + annotation_dir)
		stdout_log = os.path.join(annotation_dir, '%s.prokka.bacteria.stdout.log' % contig_base)
		score = os.path.join(annotation_dir, '%s.gene_score.txt' % contig_base)
		gene_nuc = os.path.join(annotation_dir, '%s.ffn' % contig_base)
		gene_aa = os.path.join(annotation_dir, '%s.faa' % contig_base)
		gff_file = os.path.join(annotation_dir, '%s.gff' % contig_base)
		fna_file_tmp.append(gene_nuc)

		workflow.add_task_gridable('prokka --prefix [args[0]] --addgenes --addmrna --force --metagenome '
		                           '--cpus [args[2]] '
		                           '--outdir [args[1]] [depends[0]] '
		                           '>[args[3]] 2>&1 ',
		                           depends = [contig, TrackedExecutable("prokka")],
		                           targets = [gene_nuc, gene_aa, gff_file],
		                           args = [sample, annotation_dir, threads, stdout_log],
		                           cores = threads,
		                           mem = mem_equation,
		                           time = time_equation,
								   name = contig_base + "__prokka")
	
	for myfile in fna_file_tmp:
		mym = re.search("([^\/]+)$", myfile)
		myname = mym.group(1)
		myfile_new = prokka_dir + "/" + myname
		fna_file.append(myfile_new)
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [myfile],
			targets = [myfile_new],
			name = "ln__" + myname)
		myfaa = re.sub(".ffn", ".faa", myfile)
		myfaa_new = re.sub(".ffn", ".faa", myfile_new)
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [myfaa],
			targets = [myfaa_new],
			name = "ln__" + myname)
		mygff = re.sub(".ffn", ".gff", myfile)
		mygff_new = re.sub(".ffn", ".gff", myfile_new)
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [mygff],
			targets = [mygff_new],
			name = "ln__" + myname)
	
	
	# ================================================
	# Summarize sequences
	# ================================================
	mem_equation = "50000"
	### combine gene sequences ###
	mylog = re.sub(".fna", ".log", gene_file)
	workflow.add_task('combine_gene_sequences.py -p [args[0]] -e ffn -o [targets[0]] > [args[1]] 2>&1 ',
	                  depends=utilities.add_to_list(fna_file,TrackedExecutable("combine_gene_sequences.py")),
	                  targets=[gene_file],
	                  args=[prokka_dir, mylog],
	                  name="combine_gene_sequences")

	### combine protein sequences ###
	## collect sequences
	mylog = re.sub(".faa", ".log", protein_file)
	workflow.add_task('format_protein_sequences.py -p [args[0]] -q [args[1]] -e faa -o [targets[0]] '
	                  '-m [targets[1]] >[args[2]] 2>&1 ',
	                  depends=utilities.add_to_list(faa_file,TrackedExecutable("format_protein_sequences.py")),
	                  targets=[protein_file, gene_info],
	                  args=[prokka_dir, prodigal_dir, mylog],
	                  name="format_protein_sequences")

	## sort by length and filter out short-length sequence
	mylog = re.sub(".faa", ".log", protein_sort)
	workflow.add_task('usearch -sortbylength [depends[0]] '
	                  '-fastaout [targets[0]] -minseqlength 0 >[args[0]] 2>&1 ',
	                  depends = [protein_file, TrackedExecutable("usearch")],
	                  targets = [protein_sort],
	                  args = [mylog],
					  name = "usearch__sorting")

	## extract nucleotide sequence for protein coding genes
	mylog = re.sub(".fna", ".log", gene_PC_file)
	workflow.add_task(
		'extract_protein_coding_genes.py -g [depends[0]] -p [depends[1]] -o [targets[0]] > [args[0]] 2>&1 ',
		depends = [gene_file, protein_sort, TrackedExecutable("extract_protein_coding_genes.py")],
		targets = [gene_PC_file],
		args = [mylog],
		name = "extract_protein_coding_genes")

	## extract sequences
	mylog = re.sub(".fna", ".log", complete_gene)
	workflow.add_task(
		'extract_complete_ORF_seq.py -t complete -m [depends[0]] -i [depends[1]] -o [targets[0]] >[args[0]] 2>&1',
		depends = [gene_info, gene_PC_file, TrackedExecutable("extract_complete_ORF_seq.py")],
		targets = [complete_gene],
		args = [mylog],
		name = 'extract_complete_ORF_seq')

	mylog = re.sub(".faa", ".log", complete_protein)
	workflow.add_task(
		'extract_complete_ORF_seq.py -t complete -m [depends[0]] -i [depends[1]] -o [targets[0]] >[args[0]] 2>&1',
		depends = [gene_info, protein_sort, TrackedExecutable("extract_complete_ORF_seq.py")],
		targets = [complete_protein],
		args = [mylog],
		name = 'extract_complete_ORF_seq')
	return complete_gene, complete_protein


def gene_catalog (workflow, complete_gene, complete_protein,
                  input_dir, sample_file, file_extension, threads,
                  prefix_gene_catalog, gene_catalog, gene_catalog_nuc, gene_catalog_prot,
                  mapping_dir, gene_catalog_saf, gene_catalog_count):
	"""
    This set of tasks will build gene catalogs.

    Args:
		workflow (anadama2.workflow): An instance of the workflow class.
		complete_gene: The fasta file of gene nucleotide sequences for complete ORFs.
        complete_protein: The fasta file of protein sequences for complete ORFs.
		mapping_dir: The direcory path of mapping results.
        prefix_gene_catalog: The prefix of gene catalog file.
        gene_catalog: The gene catalog file.
        gene_catalog_nuc: The fastq file of nucleotide sequences for gene catalogs.
        gene_catalog_prot: The fastq file of protein sequences for gene catalogs.
        gene_catalog_saf: The SAF gtf file for gene catalogs.
        gene_catalog_count: The count file for gene catalogs.


    Requires:
        bowtie2 (tested with 2.3.2)
        samtools (tested with 1.5)
        featureCounts (tested with Version 1.6.2)
        the nucleotide and amino acid sequences for gene catalogs
        fastq files for each sample

    Returns:
        string: file names of gene catalogs

    Example:
        from anadama2 import Workflow
        from MetaWIBELE.characterize import characterization

        # create an anadama2 workflow instance
        workflow=Workflow()

        # add quality control tasks for the fastq files
		mygene_catalog, mycounts = preprocessing_tasks.gene_catalogs (workflow, complete_gene, complete_protein,
		                                                              mapping_dir,
		                                                              prefix_gene_catalog, gene_catalog, gene_catalog_nuc, gene_catalog_prot,
		                                                              gene_catalog_saf, gene_catalog_count)

        # run the workflow
        workflow.go()
    """

	### run gene-catalog workflow
	mylog = gene_catalog_nuc + ".log"
	myclust = gene_catalog_nuc + ".clstr"
	workflow.add_task(
		'cd-hit-est -i [depends[0]] [args[0]] -o [targets[0]] >[args[1]] 2>&1 ',
		depends = [complete_gene, TrackedExecutable("cd-hit-est")],
		targets = [gene_catalog_nuc, myclust],
		args = [config.cd_hit_gene_opts, mylog],
		cores = threads,
		name = "cd-hit-est")

	mylog = gene_catalog + ".log"
	workflow.add_task(
		'extract_cluster_CD-hit.py -c [depends[0]] -o [targets[0]] >[args[0]] 2>&1 ',
		depends = [myclust, TrackedExecutable("extract_cluster_CD-hit.py")],
		targets = [gene_catalog],
		args = [mylog],
		cores = threads,
		name = "extract_cluster_CD-hit")

	mylog = gene_catalog_prot + ".log"
	workflow.add_task(
		'extract_non_redundance_seq.py -r [depends[0]] -i [depends[1]] -o [targets[0]] >[args[0]] 2>&1 ',
		depends = [gene_catalog_nuc, complete_protein, TrackedExecutable("extract_non_redundance_seq.py")],
		targets = [gene_catalog_prot],
		args = [mylog],
		cores = threads,
		name = "extract_non_redundance_seq")

	### get the abundance of gene catalog
	# run gene-abundance workflow
	mylog = gene_catalog_saf + ".log"
	workflow.add_task(
		'gene_abundance_indexRef.py -r [depends[0]] -t gene -b [args[0]] -o [targets[0]] >[args[1]] 2>&1 ',
		depends = [gene_catalog_nuc, TrackedExecutable("gene_abundance_indexRef.py")],
		targets = [gene_catalog_saf],
		args = [prefix_gene_catalog, mylog],
		cores = threads,
		name = "gene_abundance_indexRef")

	## collect sequences
	samples = []
	open_file = open(sample_file, "r")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		samples.append(info[0])
	# foreach sample
	open_file.close()

	## bowtie2 will map reads to gene categories
	flt_seqs = []
	for sample in samples:
		seq_file = "NA"
		if file_extension != "none":
			tmp = file_extension.split(",")
			for item in tmp:
				if seq_file == "NA":
					seq_file = os.path.join(input_dir, sample + '%s' % item)
				else:
					seq_file = seq_file + "," + os.path.join(input_dir, sample + '%s' % item)
		# r_seq_pair = os.path.join(qc_out_dir, '%s_adapRev_paired_2.fastq.gz' % seq_base)
		flt_seqs.append((sample, seq_file))
	# foreah sample

	## Now run bowtie2 to map reads to gene categories
	mappings = []
	mappings_tmp = []
	mem_equation = "2*12*1024 if file_size('[depends[0]]') < 10 else 4*12*1024"
	time_equation = "2*60 if file_size('[depends[0]]') < 10 else 2*2*60"
	for (sample, seq_file) in flt_seqs:
		seq_base = sample
		mydir = mapping_dir + "/" + sample
		os.system("mkdir -p " + mydir)
		sample_counts = os.path.join(mydir, seq_base + ".sort.bed")
		stdout_log = os.path.join(mydir, '%s.mapping.stdout.log' % seq_base)
		mappings_tmp.append(sample_counts)

		workflow.add_task_gridable(
			'gene_abundance.py -r [depends[0]] -u [args[0]] -t [args[1]] -s [args[2]] -w [args[3]] '
			'> [args[4]] 2>&1 ',
			depends = [gene_catalog_nuc, TrackedExecutable("gene_abundance.py")],
			targets = [sample_counts],
			args = [seq_file, threads, seq_base, mydir, stdout_log],
			cores = threads,
			mem = mem_equation,
			time = time_equation,
			name = sample + "__gene_abundance")

	for myfile in mappings_tmp:
		mym = re.search("([^\/]+)$", myfile)
		myname = mym.group(1)
		myfile_new = mapping_dir + "/" + myname
		mappings.append(myfile_new)
		workflow.add_task(
			"ln -fs [depends[0]] [targets[0]]",
			depends = [myfile],
			targets = [myfile_new],
			name = "ln__" + myname)

	# collect abundance
	mylog = gene_catalog_count + ".log"
	workflow.add_task(
		'gene_catalog_abundance.py -p [args[0]] -s sort.bed -o [targets[0]] >[args[1]] 2>&1 ',
		depends = utilities.add_to_list(mappings,TrackedExecutable("gene_catalog_abundance.py")),
		targets = [gene_catalog_count],
		args = [mapping_dir, mylog],
		name = "gene_catalog_abundance")

	return gene_catalog, gene_catalog_count


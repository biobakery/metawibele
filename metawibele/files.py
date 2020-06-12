"""
MetaWIBELE: files module
A collection of file names used by tasks

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
import copy
import sys

from anadama2 import reporters

from .utilities import name_files


class FileInfo(object):
	def __init__(self, name=None, subfolder=None, tag=None, extension=None, description=None):
		# set a list of non-path keywords
		self.non_path_keywords = ["description"]

		# concat multiple strings if present in description
		if description and isinstance(description, tuple):
			description = "\n".join(description)

		keywords = {"names": name, "subfolder": subfolder, "tag": tag, "extension": extension,
		            "description": description}
		self.keywords = {key: value for key, value in keywords.items() if value}

	def get_path_keywords(self):
		info = copy.copy(self.keywords)

		# remove the non-path keywords
		for key in self.non_path_keywords:
			try:
				del info[key]
			except KeyError:
				pass

		return info

	def __getitem__(self, key):
		""" Return the file info """
		try:
			value = self.keywords[key]
		except KeyError:
			value = ""
		return value


class Workflow(object):
	file_info = {}

	file_info["log"] = FileInfo(reporters.LOG_FILE_NAME,
	                            description="The AnADAMA2 workflow log.")

	@classmethod
	def path(cls, name, main_folder="", none_if_not_found=None, error_if_not_found=None, **keywords):
		merged_keywords = copy.copy(keywords)
		merged_keywords.update(cls.file_info[name].get_path_keywords())
		file_path = name_files(folder=main_folder, **merged_keywords)

		# if set, error if the file does not exist
		if error_if_not_found and not os.path.isfile(file_path):
			message = "\nERROR: Unable to find file: " + file_path
			desc = cls.description(name)
			if desc:
				message += "\n\nFile description:\n" + desc
			sys.exit(message)

		# if set, check if the file exists, if not return None
		if none_if_not_found and not os.path.isfile(file_path):
			file_path = None

		return file_path

	@classmethod
	def description(cls, name):
		try:
			desc = cls.file_info[name].keywords["description"]
		except (KeyError, AttributeError):
			desc = ""

		return desc

	@classmethod
	def list_file_path_description(cls, folder, input_files):
		""" List the file names and descriptions in a format to be used in an argument help description """

		desc = ""
		for required in input_files:
			desc += "\n\n".join(
				["* " + cls.path(name, folder) + " ( " + required + " )\n-- " + cls.description(name) for name in
				 input_files[required]]) + "\n"

		return desc


class Characterization(Workflow):
	""" A collection of information of folders/files created by the MetaWIBELE-characterize tasks """

	file_info = copy.copy(Workflow.file_info)

	# set the folder names annotation data workflows
	annotation_name = "characterization"
	cluster_folder_name = "cluster"
	homology_folder_name = "homology_annotation"
	non_homology_folder_name = "non_homology_annotation"
	abundance_folder_name = "abundance"
	differential_abundance_folder_name = "differential_abundance"

	# set the cluster file name
	file_info["proteinfamilies"] = FileInfo("proteinfamilies.clstr",
	                                       subfolder=os.path.join(annotation_name, cluster_folder_name),
	                                       description=("A fasta formated file.",
		                                       "Each cluster is started with >, ",
		                                       "with representative and description following. ",
		                                       "Each member of this cluster is listed next."))
	file_info["proteinfamilies_seq"] = FileInfo("proteinfamilies.rep.fasta",
	                                           subfolder=os.path.join(annotation_name, cluster_folder_name),
	                                           description=("A fasta formated file for amino acid sequences of reprecentative. "))

	# set the all feature counts file names
	file_info["feature_counts"] = FileInfo("proteinfamilies_counts.tsv",
	                                       subfolder=os.path.join(annotation_name, abundance_folder_name),
	                                       description=("A tab-delimited file with protein families as rows and samples ",
		                                       "as columns. This file includes the total feature counts (non-stratified)",
		                                       "for the features computed by MetaWIBELE."))

	# set the normed feature file names
	file_info["feature_relab"] = FileInfo("proteinfamilies_relab.tsv",
	                                      subfolder=os.path.join(annotation_name, abundance_folder_name),
	                                      description=("A tab-delimited file with features as columns and samples ",
		                                      "as rows. This file is a merged set of protein family abundances for all ",
		                                      "samples computed by MetaWIBELE. This file contains non-stratified counts ",
		                                      "of relative abundances."))

	# set the names for the rna/dna normed files
	file_info["feature_norm_ratio"] = FileInfo("proteinfamilies_rna_dna_ratio.tsv",
	                                           subfolder=os.path.join(annotation_name, abundance_folder_name),
	                                           description=("A tab-delimited file with samples as columns and genes as rows. ",
		                                           "This file includes the normalized RNA abundances as a ratio to DNA abundance. ",
		                                           "This file does not include stratified features."))

	# set the homology-based annotation file names
	file_info["uniref_mapping"] = FileInfo("proteinfamilies.uniref_stat.tsv",
	                                       subfolder=os.path.join(annotation_name, homology_folder_name),
	                                       description=("A tab-delimited file with mapping statistics as columns and protein ID as rows.",
		                                       "This file contains the homology profiles computed by DIAMOND for all proteins."))
	file_info["uniref_proteinfamilies_func"] = FileInfo("summary_UniRef90_proteinfamilies.detail.tsv",
	                                                   subfolder=os.path.join(annotation_name, homology_folder_name),
	                                                   description=("A tab-delimited file with protein families as rows and functional annotation details as columns. ",
		                                                   "This file contains the homology-based annotation of protein families using the ",
		                                                   "clusters identified by Clustering method."))
	file_info["uniref_protein_func"] = FileInfo("summary_UniRef90_proteinfamilies.ORF.detail.tsv",
	                                            subfolder=os.path.join(annotation_name, homology_folder_name),
	                                            description=("A tab-delimited file with protein as rows and functional annotation details as columns. ",
		                                            "This file contains the homology-based annotation of proteins."))
	file_info["uniref_proteinfamilies_taxa"] = FileInfo("summary_proteinfamilies_annotation.uniref90_annotation.tsv",
	                                                   subfolder=os.path.join(annotation_name, homology_folder_name),
	                                                   description=("A tab-delimited file with protein families as rows and taxonomy annotations as columns. ",
		                                                   "This file contains the homology-based annotation of protein families."))
	file_info["uniref_protein_taxa"] = FileInfo("summary_protein_annotation.uniref90_annotation.tsv",
	                                            subfolder=os.path.join(annotation_name, homology_folder_name),
	                                            description=("A tab-delimited file with protein as rows and taxonomy annotations as columns. ",
		                                            "This file contains the homology-based annotation of proteins."))
	file_info["BGC_proteinfamilies"] = FileInfo("summary_antiSMASH_proteinfamilies.detail.tsv",
	                                           subfolder=os.path.join(annotation_name, homology_folder_name,
	                                                                  "InterProScan"),
	                                           description=("A tab-delimited file with protein families as rows and antiSMASH annotation details as columns. ",
		                                           "This file contains the antiSMASH annotation of protein families using the ",
		                                           "clusters identified by Clustering method."))
	file_info["BGC_protein"] = FileInfo("summary_antiSMASH_proteinfamilies.ORF.detail.tsv",
	                                    subfolder=os.path.join(annotation_name, homology_folder_name,
	                                                           "InterProScan"),
	                                    description=("A tab-delimited file with protein as rows and antiSMASH annotation details as columns. ",
		                                    "This file contains the antiSMASH annotation of proteins"))

	# set the non-homology-based annotation file names
	file_info["signalp_proteinfamilies"] = FileInfo("summary_SignalP_proteinfamilies.detail.tsv",
	                                               subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                      "InterProScan"),
	                                               description=("A tab-delimited file with protein families as rows and signalp annotation details as columns. ",
		                                               "This file contains the signalp annotation of protein families using the ",
		                                               "clusters identified by Clustering method."))
	file_info["signalp_protein"] = FileInfo("summary_SignalP_proteinfamilies.ORF.detail.tsv",
	                                        subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                               "InterProScan"),
	                                        description=("A tab-delimited file with protein as rows and signalp annotation details as columns. ",
		                                        "This file contains the signalp annotation of proteins."))
	file_info["tmhmm_proteinfamilies"] = FileInfo("summary_TMHMM_proteinfamilies.detail.tsv",
	                                             subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                    "InterProScan"),
	                                             description=("A tab-delimited file with protein families as rows and tmhmm annotation details as columns. ",
		                                             "This file contains the tmhmm annotation of protein families using the ",
		                                             "clusters identified by Clustering method."))
	file_info["tmhmm_protein"] = FileInfo("summary_TMHMM_proteinfamilies.ORF.detail.tsv",
	                                      subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                             "InterProScan"),
	                                      description=("A tab-delimited file with protein as rows and tmhmm annotation details as columns. ",
		                                      "This file contains the tmhmm annotation of proteins."))
	file_info["phobius_proteinfamilies"] = FileInfo("summary_Phobius_proteinfamilies.detail.tsv",
	                                               subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                      "InterProScan"),
	                                               description=("A tab-delimited file with protein families as rows and Phobius annotation details as columns. ",
		                                               "This file contains the Phobius annotation of protein families using the ",
		                                               "clusters identified by Clustering method."))
	file_info["phobius_protein"] = FileInfo("summary_Phobius_proteinfamilies.ORF.detail.tsv",
	                                        subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                               "InterProScan"),
	                                        description=("A tab-delimited file with protein as rows and Phobius annotation details as columns. ",
		                                        "This file contains the Phobius annotation of proteins."))
	file_info["denovo_proteinfamilies_signal"] = FileInfo("summary_Denovo_proteinfamilies.signaling.detail.tsv",
	                                                     subfolder=os.path.join(annotation_name,
	                                                                            non_homology_folder_name,
	                                                                            "InterProScan"),
	                                                     description=("A tab-delimited file with protein families as rows and signaling annotation details as columns. ",
		                                                     "This file contains the signaling annotation of protein families using the ",
		                                                     "clusters identified by Clustering method."))
	file_info["denovo_protein_signal"] = FileInfo("summary_Denovo_proteinfamilies.ORF.signaling.detail.tsv",
	                                              subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                     "InterProScan"),
	                                              description=("A tab-delimited file with protein as rows and signaling annotation details as columns. ",
		                                              "This file contains the signaling annotation of proteins."))
	file_info["denovo_proteinfamilies_transmembrane"] = FileInfo("summary_Denovo_proteinfamilies.transmembrane.detail.tsv",
																subfolder=os.path.join(annotation_name, non_homology_folder_name,
		                                                                                "InterProScan"),
																description=("A tab-delimited file with protein families as rows and transmembrane annotation details as columns. ",
																	"This file contains the transmembrane annotation of protein families using the ",
																	"clusters identified by Clustering method."))
	file_info["denovo_protein_transmembrane"] = FileInfo("summary_Denovo_proteinfamilies.ORF.transmembrane.detail.tsv",
	                                                     subfolder=os.path.join(annotation_name,
	                                                                            non_homology_folder_name,
	                                                                            "InterProScan"),
	                                                     description=("A tab-delimited file with protein as rows and transmembrane annotation details as columns. ",
		                                                     "This file contains the transmembrane annotation of proteins."))
	file_info["interpro_proteinfamilies"] = FileInfo("summary_InterProScan_proteinfamilies.detail.tsv",
	                                                subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                       "InterProScan"),
	                                                description=("A tab-delimited file with protein families as rows and InterProScan annotation details as columns. ",
		                                                "This file contains the InterProScan annotation of protein families using the ",
		                                                "clusters identified by Clustering method."))
	file_info["interpro_protein"] = FileInfo("summary_InterProScan_proteinfamilies.ORF.detail.tsv",
	                                         subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                "InterProScan"),
	                                         description=("A tab-delimited file with protein as rows and InterProScan annotation details as columns. ",
		                                         "This file contains the InterProScan annotation of proteins."))
	file_info["pfam_proteinfamilies"] = FileInfo("summary_Pfam_proteinfamilies.detail.tsv",
	                                            subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                   "Pfam"),
	                                            description=("A tab-delimited file with protein families as rows and Pfam annotation details as columns. ",
		                                            "This file contains the Pfam annotation of protein families using the ",
		                                            "clusters identified by Clustering method."))
	file_info["pfam_protein"] = FileInfo("summary_Pfam_proteinfamilies.ORF.detail.tsv",
	                                     subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                            "Pfam"),
	                                     description=("A tab-delimited file with protein as rows and Pfam annotation details as columns. ",
		                                     "This file contains the Pfam annotation of proteins."))
	file_info["pfam2go_proteinfamilies"] = FileInfo("summary_Pfam2GO_proteinfamilies.detail.tsv",
	                                               subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                      "Pfam"),
	                                               description=("A tab-delimited file with protein families as rows and Pfam2GO annotation details as columns. ",
		                                               "This file contains the Pfam2GO annotation of protein families using the ",
		                                               "clusters identified by Clustering method."))
	file_info["pfam2go_protein"] = FileInfo("summary_Pfam2GO_proteinfamilies.ORF.detail.tsv",
	                                        subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                               "Pfam"),
	                                        description=("A tab-delimited file with protein as rows and Pfam2GO annotation details as columns. ",
		                                        "This file contains the Pfam2GO annotation of proteins."))
	file_info["DDI_proteinfamilies"] = FileInfo("summary_DOMINE_proteinfamilies.detail.tsv",
	                                           subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                  "Pfam"),
	                                           description=("A tab-delimited file with protein families as rows and DOMINE annotation details as columns. ",
		                                           "This file contains the DOMINE annotation of protein families using the ",
		                                           "clusters identified by Clustering method."))
	file_info["DDI_protein"] = FileInfo("summary_DOMINE_proteinfamilies.ORF.detail.tsv",
	                                    subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                           "Pfam"),
	                                    description=("A tab-delimited file with protein as rows and DOMINE annotation details as columns. ",
		                                    "This file contains the DOMINE annotation of proteins."))
	file_info["SIFTS_proteinfamilies"] = FileInfo("summary_SIFTS_proteinfamilies.detail.tsv",
	                                             subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                    "Pfam"),
	                                             description=("A tab-delimited file with protein families as rows and SIFTS annotation details as columns. ",
		                                             "This file contains the SIFTS annotation of protein families using the ",
		                                             "clusters identified by Clustering method."))
	file_info["SIFTS_protein"] = FileInfo("summary_SIFTS_proteinfamilies.ORF.detail.tsv",
	                                      subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                             "Pfam"),
	                                      description=("A tab-delimited file with protein as rows and SIFTS annotation details as columns. ",
		                                      "This file contains the SIFTS annotation of proteins."))
	file_info["ExpAtlas_proteinfamilies"] = FileInfo("summary_ExpAtlas_proteinfamilies.detail.tsv",
	                                                subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                       "Pfam"),
	                                                description=("A tab-delimited file with protein families as rows and ExpAtlas annotation details as columns. ",
		                                                "This file contains the ExpAtlas annotation of protein families using the ",
		                                                "clusters identified by Clustering method."))
	file_info["ExpAtlas_protein"] = FileInfo("summary_ExpAtlas_proteinfamilies.ORF.detail.tsv",
	                                         subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                "Pfam"),
	                                         description=("A tab-delimited file with protein as rows and ExpAtlas annotation details as columns. ",
		                                         "This file contains the ExpAtlas annotation of proteins."))
	file_info["psortb_proteinfamilies"] = FileInfo("summary_PSORTb_proteinfamilies.detail.tsv",
	                                              subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                     "PSORTb"),
	                                              description=("A tab-delimited file with protein families as rows and PSORTb annotation details as columns. ",
		                                              "This file contains the PSORTb annotation of protein families using the ",
		                                              "clusters identified by Clustering method."))
	file_info["psortb_protein"] = FileInfo("summary_PSORTb_proteinfamilies.ORF.detail.tsv",
	                                       subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                              "PSORTb"),
	                                       description=("A tab-delimited file with protein as rows and PSORTb annotation details as columns. ",
		                                       "This file contains the PSORTb annotation of proteins."))
	file_info["dna_abu_proteinfamilies"] = FileInfo("summary_DNA_proteinfamilies.abundance.detail.tsv",
	                                               subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                      abundance_folder_name),
	                                               description=("A tab-delimited file with protein families as rows and DNA-abundance annotation details as columns. ",
		                                               "This file contains the DNA-abundance annotation of protein families using the ",
		                                               "clusters identified by Clustering method."))
	file_info["dna_abu_protein"] = FileInfo("summary_DNA_proteinfamilies.ORF.abundance.detail.tsv",
	                                        subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                               abundance_folder_name),
	                                        description=("A tab-delimited file with protein as rows and DNA-abundance annotation details as columns. ",
		                                        "This file contains the DNA-abundance annotation of proteins."))
	file_info["rna_abu_proteinfamilies"] = FileInfo("summary_RNA_proteinfamilies.abundance.detail.tsv",
	                                               subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                      abundance_folder_name),
	                                               description=("A tab-delimited file with protein families as rows and RNA-abundance annotation details as columns. ",
		                                               "This file contains the RNA-abundance annotation of protein families using the ",
		                                               "clusters identified by Clustering method."))
	file_info["rna_abu_protein"] = FileInfo("summary_RNA_proteinfamilies.ORF.abundance.detail.tsv",
	                                        subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                               abundance_folder_name),
	                                        description=("A tab-delimited file with protein as rows and RNA-abundance annotation details as columns. ",
		                                        "This file contains the RNA-abundance annotation of proteins."))
	file_info["msp_proteinfamilies"] = FileInfo("summary_MSPminer_proteinfamilies.abundance.detail.tsv",
	                                           subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                                  abundance_folder_name),
	                                           description=("A tab-delimited file with protein families as rows and MSP annotation details as columns. ",
		                                           "This file contains the MSP annotation of protein families using the ",
		                                           "clusters identified by Clustering method."))
	file_info["msp_protein"] = FileInfo("summary_MSPminer_proteinfamilies.ORF.abundance.detail.tsv",
	                                    subfolder=os.path.join(annotation_name, non_homology_folder_name,
	                                                           abundance_folder_name),
	                                    description=("A tab-delimited file with protein as rows and MSP annotation details as columns. ",
		                                    "This file contains the MSP annotation of proteins."))
	file_info["proteinfamilies_taxa"] = FileInfo("summary_proteinfamilies_annotation.taxonomy.tsv",
	                                                   subfolder=os.path.join(annotation_name, homology_folder_name),
	                                                   description=("A tab-delimited file with protein families as rows and taxonomy annotations as columns. ",
		                                                   "This file contains the non-homology based annotation of protein families."))
	file_info["protein_taxa"] = FileInfo("summary_protein_annotation.taxonomy.tsv",
	                                            subfolder=os.path.join(annotation_name, homology_folder_name),
	                                            description=("A tab-delimited file with protein as rows and taxonomy annotations as columns. ",
		                                            "This file contains the non-homology based annotation of proteins."))

	# set the differential abundance file names
	file_info["dna_DA_proteinfamilies"] = FileInfo("summary_DNA_proteinfamilies.DA.detail.tsv",
												subfolder=os.path.join(annotation_name, differential_abundance_folder_name,
																		"DNA"),
												description = ("A tab-delimited file with protein families as rows and DNA-based differential abundance annotation details as columns. ",
													"This file contains the DNA-based differential abundance annotation of protein families using the ",
													"clusters identified by Clustering method."))
	file_info["dna_DA_protein"] = FileInfo("summary_DNA_proteinfamilies.ORF.DA.detail.tsv",
											subfolder=os.path.join(annotation_name, differential_abundance_folder_name,
											                       "DNA"),
											description = ("A tab-delimited file with protein as rows and DNA-based differential abundance annotation details as columns. ",
												"This file contains the DNA-based differential abundance annotation of proteins."))
	file_info["rna_DE_proteinfamilies"] = FileInfo("summary_RNA_proteinfamilies.DE.detail.tsv",
												subfolder=os.path.join(annotation_name, differential_abundance_folder_name,
	                                                                     "RNA"),
												description = ("A tab-delimited file with protein families as rows and RNA-based differential expression annotation details as columns. ",
													"This file contains the DNA-based differential expression annotation of protein families using the ",
													"clusters identified by Clustering method."))
	file_info["rna_DE_protein"] = FileInfo("summary_RNA_proteinfamilies.ORF.DE.detail.tsv",
											subfolder=os.path.join(annotation_name, differential_abundance_folder_name,
	                                                              "RNA"),
											description = ("A tab-delimited file with protein as rows and RNA-based differential expression annotation details as columns. ",
												"This file contains the RNA-based differential expression annotation of proteins"))

	# set name for finalized characterization files
	file_info["annotation_table"] = FileInfo("proteinfamilies_annotation_table.tsv",
	                                        subfolder=os.path.join(annotation_name),
											description = ("A tab-delimited file with protein families as rows and formatted annotation as columns. ",
												"This file contains the annotation of protein families using the ",
												"clusters identified by Clustering method."))
	file_info["annotation_table"] = FileInfo("attribute_table.tsv",
	                                        subfolder=os.path.join(annotation_name),
											description = ("A tab-delimited file with protein families as rows and attributes as columns. ",
												"This file contains the specific annotations of protein families using differen methods."))


	class Prioritization(Workflow):
		""" A collection of information of folders/files created by the MetaWIBELE-prioritize tasks"""

		file_info = copy.copy(Workflow.file_info)

		# set the folder names prioritization data workflows
		annotation_name = "prioritization"

		file_info["unsupervised_score"] = FileInfo("unsupervised_prioritization_score_table.tsv",
		                                           description="A tab-delimited file with protein families as rows and ranking as " + \
		                                                       "columns. Includes the prioritization using abundance and prevalence profiles.")
		file_info["supervised_score"] = FileInfo("supervised_prioritization_score_table.tsv",
		                                         description="A tab-delimited file with protein families as rows and ranking as " + \
		                                                     "columns. Includes the prioritization using abundance, prevalence profiles " + \
		                                                     "and association with phenotypes.")

		file_info["unsupervised_table"] = FileInfo("unsupervised_prioritization_table.tsv",
		                                           description="A tab-delimited file with protein families as rows and ranking as " + \
		                                                       "columns. Includes the prioritization using abundance and prevalence profiles.")
		file_info["supervised_table"] = FileInfo("supervised_prioritization_table.tsv",
		                                         description="A tab-delimited file with protein families as rows and ranking as " + \
		                                                     "columns. Includes the prioritization using abundance, prevalence profiles " + \
		                                                     "and association with phenotypes.")

		file_info["interested_table"] = FileInfo("interested_prioritization_table.tsv",
		                                         description="A tab-delimited file with protein families as rows and ranking as " + \
		                                                     "columns. Includes the prioritization using abundance, prevalence profiles " + \
		                                                     "and association with phenotypes.")

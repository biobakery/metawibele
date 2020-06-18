# MetaWIBELE User Manual

**MetaWIBELE** (**W**orkflow to **I**dentify novel **B**ioactive **Ele**ments in the microbiome) is a computational pipeline that identifies novel bioactive microbial gene products from metagenomes and finds new immunomodulatory gene families, especially targeting secreted/extracellular proteins to enrich for likely host interactors. The prioritized list of gene products can be further used for downstream experimental validation. MetaWIBELE is available as module of bioBakery [bioBakery repository](https://github.com/biobakery).

## Citing MetaWIBELE

**A manuscript describing MetaWIBELE is currently in prep:**

Identifying Novel Bioactive Microbial Gene Products in Inflammatory Bowel Disease

**And feel free to link to MetaWIBELE in your Methods:**

[http://huttenhower.sph.harvard.edu/MetaWIBELE](http://huttenhower.sph.harvard.edu/MetaWIBELE)

**For additional information, read the** [MetaWIBELE Tutorial](https://github.com/biobakery/metawibele)

Support for MetaWIBELE is available via [the MetaWIBELE channel](https://forum.biobakery.org/c/Microbial-community-profiling/MetaWIBELE) of the bioBakery Support Forum.

***

## Contents ##

* [Workflow](#workflow)
	* [Workflow by bypass mode](#workflow-by-bypass-mode) 
* [Install MetaWIBELE](#install-metawibele)
    * [Requirements](#requirements)
    * [Installation](#installation)
    	* [Download MetaWIBELE](#download-metawibele)
    	* [Install MetaWIBELE](#install-metawibele)
    	* [Download databases](#download-databases)
    	* [Download configuration files](#download-configuration-files)
			* [Download global configuration template](#download-global-configuration-template)
    		* [Download local configuration template](#download-local-configuration-template)
    		* [Download vignette configuraion template](#download-vignette-configuration-template)
* [Quick-start Guide](#quick-start-guide)
    * [How to run](#how-to-run)
    * [Standard Workflows](#standard-workflows)
    	* [MetaWIBELE-characterize workflow](#metawibele-characterize-workflow)
    		* [Input files for characterization](#input-files-for-characterization)
    		* [MetaWIBELE-characterize workflow](#metawibele-characterize-workflow)
    		* [Demo run of MetaWIBELE-characterize](#demo-run-of-metawibele-characterize)
    		* [Output files of MetaWIBELE-characterize](#output-files-of-metawibele-characterize)
    	* [MetaWIBELE-prioritize workflow](#metawibele-prioritize-workflow)
    		* [Input files for prioritization](#input-files-for-prioritization)
    		* [MetaWIBELE-prioritize workflow](#metawibele-prioritize-workflow)
    		* [Demo run of MetaWIBELE-prioritize](#demo-run-of-metawibele-prioritize)
    		* [Output files of MetaWIBELE-prioritize](#output-files-of-metawibele-prioritize)
* [Guides to MetaWIBELE Utilities](#guides-to-metawibele-utilities)
	* [Preprocessing sequencing reads to build gene catalogs](#preprocessing-sequencing-reads-into-to-build-gene-catalogs)
		* [Preprocessing workflow](#preprocessing-workflow)
		* [Input files for preprocessing workflow](#input-files-preprocessing-workflow)
		* [Demo run of preprocessing workflow](#demo-run-of-preprocessing-workflow) 
		* [Output files of preprocessing workflow](#output-files-of-preprocessing-workflow)
* [Download MetaWIBELE resources](#download-metawibele-resources)
	* [Information of gene catalogs](#information-of-gene-catalogs)
	* [Characterization of protein families](#characterization-of-protein-families)
	* [Prioritization of protein families](#prioritization-of-protein-families)
    	
***


## Workflow
![workflow.png](https://www.dropbox.com/s/vsg7baww6utske1/MetaWIBELE_overview.png?raw=1)
***

### Workflow by bypass mode
There are multiple bypass options that will allow you to adjust the standard workflow.

Bypass options:

* --bypass-global-homology
	* do not annotate protein families based on global homology information
* --bypass-domain-motif
	* do not annotate protein families based on domain/motif information
* --bypass-abundance
	* do not annotate protein families based on abundance information
* --bypass-integration
	* do not integrate annotations for protein families
* --bypass-optional
	* do not prioritize protein families based on selecting our for interested annotations (optional prioritization)


## Install MetaWIBELE
### Requirements

1. [Python](https://www.python.org/) (version >= 3.7)
2. [AnADAMA2](https://huttenhower.sph.harvard.edu/anadama2) (version >= 0.7.4-devel)
3. [CD-hit](http://weizhongli-lab.org/cd-hit/) (version >= 4.7)
4. [Diamond](http://www.diamondsearch.org/index.php) (version >= 0.9.5)
5. [MSPminer](https://www.enterome.com/downloads/) (version >= 2)
6. [MaAsLin2](https://huttenhower.sph.harvard.edu/maaslin2) (versuib >= 1.1.2) (only required if using MaAsLin2 to associate with host phenotypes)
7. [Interproscan](https://github.com/ebi-pf-team/interproscan/wiki) (version >= 5.31-70) (only required if using Interproscan to annotate domains and motifs)
8. [Signalp](http://www.cbs.dtu.dk/services/SignalP-4.1/) (version >= 4.1) (only required if using Signalp to annotate signal peptides)
9. [TMHMM](http://www.cbs.dtu.dk/services/TMHMM/) (version >= 2.0) (only required if using TMHMM to annotate transmembrane proteins)
10. [Phobius](http://phobius.sbc.su.se/) (version >= 1.01) (only required if using Phobius to annotate both signal peptides and transmembrane proteins)
11. [PSORTb](https://psort.org/documentation/index.html) (version >= 3.0) (only required if using PSORTb to predict subcellular localization)
12. **Optional**: only required if using MetaWIBELE utitlity to preprocess metagenomic sequencing reads
	* [MEGAHIT](https://github.com/voutcn/megahit) (version >= 1.1.3) 
	* [Prokka](https://github.com/tseemann/prokka) (version >= 1.14-dev; recommend to close '-c' parameter in setting prodigal parameters)
	* [Prodigal](https://github.com/hyattpd/Prodigal) (version >= 2.6)
	* [USEARCH](http://www.drive5.com/usearch/) (version >= 9.0.2132_i86linux64)
	* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.3.2)
	* [SAMtools](https://github.com/samtools/) (version >= 1.5)
	* [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) (version >= 1.6.2)

**Note:** Please install the required software in a location in your `$PATH`. If you always run with gene catalogs, the optional softwares are not required. Also if you always run with one or more bypass options (for information on bypass options, see optional arguments to the section [Workflow by bypass mode](#workflow-by-bypass-mode)), the software required for the steps you bypass does not need to be installed.


### Installation
#### Download MetaWIBELE
You can download the latest MetaWIBELE release or the development version. The source contains example files.

Option 1: Latest Release (Recommended)

* download [MetaWIBELE.tar.gz](http://huttenhower.sph.harvard.edu/MetaWIBELE/MetaWIBELE.tar.gz) and unpack the latested release of MetaWIBELE

Option 2: Development Version

* Create a clone of the Git repository 
	* `$ git clone https://github.com/biobakery/metawibele.git`
* You can always update to the latest version of the repository with:
	* `$ git pull --update`

#### Install MetaWIBELE

* Installing from source
	* Move to the MetaWIBELE directory
		* `$ cd $MetaWIBELE_PATH`

	* Install MetaWIBELE package
		* `$ python setup.py install`
		* If you do not have write permissions to `/usr/lib/`, then add the option --user to the install command. This will install the python package into subdirectories of `~/.local/`. Please note when using the --user install option on some platforms, you might need to add `~/.local/bin/` to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `metawibele_workflow: command not found` when trying to run MetaWIBELE after installing with the --user option. 
		* Similarly, you can also specify the installation directory using --prefix option which will install the python package into the directory `$YOUR_INSTALL_DIR` that is a directory on PYTHONPATH and which Python reads ".pth" files from. You might need to add `$YOUR_INSTALL_DIR/bin` to your `$PATH` as it might not be included by default.


#### Download databases
UniRef databases are required if you will use MetaWIBELE to do global-homology based annotation and taxonomic annotation. These databases have been host by [HUMAnN](https://github.com/biobakery/humann). You can download these databases by using HUMAnN utility scripts and provid `$UNIREF_LOCATION` as the location to install the database.

* Download the full UniRef90 database (11.0GB, recommend):

	`$ humann2_databases --download uniref uniref90_diamond $UNIREF_LOCATION`

* Download additional UniRef90 mapping files (xxx GB):

	`$ humann2_databases --download utility_mapping full $UNIREF_LOCATION`

	* Mappings are available for the UniRef90 gene families to the following systems:
		* MetaCyc Reactions
		* KEGG Orthogroups (KOs)
		* Pfam domains
		* Level-4 enzyme commission (EC) categories
		* EggNOG (including COGs)
		* Gene Ontology (GO)
		* Taxonomic information of the latest common acestor (LCA)
		* Taxonomic information of the representative
		* UniProt ID corresponding to the UniRef representative
		* Protein names of the UniRef representative
		* Gene names of the UniRef representative 
	* In most cases, mappings are directly inferred from the annotation of the corresponding UniRef centroid sequence in UniProt.


	
### Download configuration files

#### Download global configuration template

To run MetaWIBELE, one global configuation file `metawibele.cfg` is required to make basic settings. 

* Download `metawibele.cfg` into  your working directory:
	* `$ metawibele_download_config --config-type global`
* Make your own configurations:
	* Specify your input files and output folder:

	```
	[input]
	# Study name
	study = 
	# Metadata file
	metadata = 
	# Sample list file
	sample_list = 
	# The clustering information of non-reduandant gene catalogs saved with an extended fasta file, similar with CD-hit output
	gene_catalog = 
	# The protein sequences of representatives of gene catalogs
	gene_catalog_prot = 
	# The reads counts matrix table across samples of gene catalogs
	gene_catalog_count = 

	[output]
	# The prefix name for output results
	basename = 
	# The output directory
	output_dir =
	``` 

	* Specify applied computational resources:

	```
	[computation]
	# The number of cores that you’re requesting [ Default: 1 ]
	threads = 1
	# The amount of memory (in MB) that you will be using for your job [ Default: 20000 ] 
	memory = 20000
	# The amount of time (in minute) that you will be using for your job [ Default: 60 ]
	time = 60
	```
	
	* Specify the path of uniref databases:
	
	```
	[uniref]
	# The uniref databases used by MetaWIBELE. [data_path] provide the absolute path of the uniref databases folder
	database = 
	```

	* Specify parameter setteings for annotations

	```
	[abundance]
	# The config file used by MSPminer. [config_file] provide the mspminer config file; [none] use the dedault config files installed in the metawibele package [ Default: none ]
	mspminer = none
	# The method for normalization [Choices: cpm, relab]. [cpm] copies per million units (sum to 1 million); [relab] relative abundance (sum to 1) [ Default: cpm ]  
	normalize = cpm
	# The minimum abundance for each feature [ Default: 0 ]   
	abundance_detection_level = 0

	[interproscan]
	# The command of interproscan, e.g. /my/path/interproscan/interproscan.sh
	interproscan_cmmd = 
	# The appls used by interroiscan: [appls] comma separated list of analyses, [ Choices: CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM,Phobius,SignalP,TMHMM ]; [none] use all all analyses for running [ Default: Pfam,Phobius,SignalP,TMHMM ]
	interproscan_appl = "Pfam,Phobius,SignalP,TMHMM"
	# The number of spliting files which can be annotated in parallel 	[ Default: 1 ]
	split_number = 1
	
	[maaslin2]
	# The command of Maaslin2, e.g. /my/path/Maaslin2/R/Maaslin2.R
	maaslin2_cmmd =
	# The minimum abundance for each feature [ Default: 0 ]  
	min_abundance = 0 
	# The minimum percent of samples for which a feature is detected at minimum abundance [ Default: 0.1 ]
	min_prevalence = 0.1  
	# The q-value threshold for significance [ Default: 0.25 ]
	max_significance = 0.25
	# The normalization method to apply [ Choices: TSS, CLR, CSS, NONE, TMM ], [ Default: TSS ]
	normalization = NONE
	# The transform to apply [ Choices: LOG, LOGIT, AST, NONE ],  [ Default: LOG ]
	transform = LOG
	# The analysis method to apply [ Choices: LM, CPLM, ZICP, NEGBIN, ZINB ], [ Default: LM ]
	analysis_method = LM
	# The fixed effects for the model, comma-delimited for multiple effects [ Default: all ]
	fixed_effects = all
	# The random effects for the model, comma-delimited for multiple effects [ Default: none ]
	random_effects = none
	# The correction method for computing the q-value [ Default: BH ]
	correction = BH
	# Apply z-score so continuous metadata are on the same scale [ Default: TRUE ]pply z-score so continuous metadata are on the same scale [ Default: TRUE ]
	standardize = TRUE
	# Generate a heatmap for the significant associations [ Default: FALSE ]
	plot_heatmap = FALSE
	# In heatmap, plot top N features with significant associations [ Default: FALSE ]
	heatmap_first_n = FALSE
	# Generate scatter plots for the significant associations [ Default: FALSE ]
	plot_scatter = FALSE
	# The number of R processes to run in parallel [ Default: 1 ]
	maaslin2_cores = 1
	# The minimum percent of case-control samples used for comparision in which a feature is detected [ Default: 0.1 ]
	tshld_prevalence = 0.10
	# The q-value threshold for significance used as DA annotations [ Default: 0.05 ]
	tshld_qvalue = 0.05
	# The statistic used as effect size [ Choices: coef, mean(log) ]. [coef] represents the coefficient from the model; [mean(log)] represents the difference of mean values between case and control conditions. [  Default: mean(log) ]
	effect_size = mean(log)
	# The main phenotype metadata used for prioritization, e.g. metadata1
	phenotype =
	# Set metadata value used as control condition for phenotype metadata varibles; use semicolon to seperate variables, e.g. "metadata1:control_status1;metadata2:control_status2"
	flag_ref = 
	# Case and control metadata pairs for phenotype metadata variables; use semicolon to seperate variables, e.g. "metadata1:case_status1|control_status1;metadata2:case_status2|control_status2,case_status3|control_status2"
	case_control_status =
	```


#### Download local configuration template
By default, MetaWIBELE will perform by using the local configuration files installed in the package. Optionally, you can also make your own local configuration files and provide them with optional arguments to MetaWIBELE. For example, the local characterization configuration file can be provided with `--characterization-config characterization.cfg`

* Download local configuration template files into your working directory:
	* `$ metawibele_download_config --config-type local`
* Mofidy and provide your own configurations:
	* Configurations for characterization in `characterization.cfg` which can be provided with `--characterization-config characterization.cfg`:
	
	```
	[global_homology]
	# protein family annotation based on global similarity: [yes] process this step, [no] skip this step [ Default: yes ]
	uniref = yes

	[domain_motif]
	# domain annotation: [yes] process this step, [no] skip this step [ Default: yes ]
	interproscan = yes
	# Pfam2GO to annotate GOs: [yes] process this step, [no] skip this step [ Default: yes ]
	pfam2go = yes
	# domain-domain interaction from DOMINE database: [yes] process this step, [no] skip this step [ Default: yes ]
	domine = yes
	# DDI with SIFTS evidence: [yes] process this step, [no] skip this step [ Default: yes ]
	sifts = yes
	# DDI with human expression from ExpAtlas database: [yes] process this step, [no] skip this step [ Default: yes ]
	expatlas = yes
	# subcellular annotation: [yes] process this step, [no] skip this step [ Default: yes ]
	psortb = yes

	[abundance]
	# summary DNA abundance: [label] provide label for DNA abundance, e.g. DNA_abundance, [no] skip this step [ Default: DNA_abundance ]
	dna_abundance = DNA_abundance
	# differential abundance based on DNA abundance: [label] provide label for DA annotation, e.g. MaAsLin2_DA, [no] skip this step [ Default: MaAsLin2_DA ]
	dna_da = MaAsLin2_DA
	
	[integration]
	# summarize annotation info: [yes] process this step, [no] skip this step [ Default: yes ]
	summary_ann = yes
	# generate finalized annotations: [yes] process this step, [no] skip this step [ Default: yes ]
	finalization = yes
	``` 
	
	* Configurations for prioritization in `prioritization.cfg` which can be provided with `--prioritization-config prioritization.cfg`:
	
	```
	## Mandatory ranking
	[unsupervised]
	# Weight value of prevalence to caculate weighted harmonic mean, named as beta parameter[ Default: 0.50 ] 
	DNA_prevalence = 0.50
	# Weight value of mean abundance to calculate weighted harmonic mean [ Default: 0.50 ] 
	DNA_abundance = 0.50

	[supervised]
	# Use the ecological property (abundance) to do prioritization. [required] required item, [optional] optional item, [none] ignoring. [ Default: required]
	DNA-within-phenotype_abundance = required
	# Use the ecological property (prevalence) to do prioritization. [required] required item, [optional] optional item, [none] ignoring. [ Default: required]
	DNA-within-phenotype_prevalence = required
	# Use the association with phenotypes (q values from associations) to do prioritization. [required] required item, [optional] optional item, [none] ignoring. [ Default: required]
	MaAsLin2_DA__qvalue = required
	# Use the association with phenotypes (effect size from associations) to do prioritization. [required] required item, [optional] optional item, [none] ignoring. [ Default: required]
	MaAsLin2_DA__mean(log) = required
	
	## Binary filtering for selection subset
	# All [vignette_type] should be true
	# All [required] items should be true 
	# At least one [optional] item should be true 
	# All [none] items will be ignored
	# Default: select protein families significantly associated with the main clinical phenotype

	[filtering]
	# Filter for interested functional vignettes type [Choices: pilin | secreted_system | other user defined | none]
	vignettes = none

	# Filter for significant associations: [required] required item, [optional] optional item, [none] ignoring [ Default: required ]
	MaAsLin2_DA-sig = required
	
	# Filter for biochemical annotations: [required] required item, [optional] optional item, [none] ignoring
	ExpAtlas_interaction = optional
	DOMINE_interaction = optional
	SIFTS_interaction = optional
	Denovo_signaling = optional
	Denovo_transmembrane = optional
	PSORTb_extracellular = optional
	PSORTb_cellWall = optional
	PSORTb_outerMembrane = optional
	UniRef90_extracellular = optional
	UniRef90_signaling = optional
	UniRef90_transmembrane = optional
	UniRef90_cellWall = optional
	UniRef90_outerMembrane = optional
	UniRef90_PfamDomain = optional
	InterProScan_PfamDomain = optional
	InterProScan_SUPERFAMILY = optional
	InterProScan_ProSiteProfiles = optional 
	InterProScan_ProSitePatterns = optional
	InterProScan_Gene3D = optional
	InterProScan_PANTHER = optional
	InterProScan_TIGRFAM = optional
	InterProScan_SFLD = optional
	InterProScan_ProDom = optional
	InterProScan_Hamap = optional
	InterProScan_SMART = optional
	InterProScan_CDD = optional
	InterProScan_PRINTS = optional
	InterProScan_PIRSF = optional
	InterProScan_MobiDBLite = optional

	```
	
#### Download vignette configuration template
MetaWIBELE can accept user defined vignette functions of interest for further prioritization. You can make your own vignettes configuration files and provide them with optional argument to MetaWIBELE. For example, the vignette function configurations file can be provided with `--vignette-config vignettes_function.tsv`

* Download local vignettes template file into your working directory:
	* `$ metawibele_download_config --config-type vignette`
* Make your own configurations:
	* This is a tab-separated values file.
	* Two required columns: `type` indicates which type of function it is; `annotation` indicates the specific annotations assigned by MetaWIBELE given a annotation type.
	* Other optional columns: `annotation_type` indicates what tye of annotation it is in MetaWIBELE; `description` indicates detailed descriptions of the annotation.
	
	```
	type    annotation       annotation_type  description
	pilin   PF11530 PfamDomain      Minor type IV pilin, PilX-like
	pilin   PF14245 PfamDomain      Type IV pilin PilA
	pilin   PF16734 PfamDomain      Type IV pilin-like putative secretion pathway protein G/H
	pilin   PF08805 PfamDomain      Type 4 secretion system, PilS, N-terminal
	pilin   PF09160 PfamDomain      FimH, mannose-binding domain
	```


## Quick-start Guide
### How to run
* For a list of all available workflows, run:

	`$ metawibele_workflow --help`
	
* All workflows follow the general command format:

	`$ metawibele_workflow $WORKFLOW --input $INPUT --output $OUTPUT`

* For specific options for a workflow, run:

	`$ metawibele_workflow $WORKFLOW --help`

* Run **characterization** workflow

	```
	metawibele_workflow characterize \
		--characterization-config $characterization_conf \
		--input $INPUT \
		--output $OUTPUT
	```

* Run **prioritization** workflow

	```
	metawibele_workflow prioritize \
		--prioritization-config $prioritization_conf \
 		--input $INPUT \
 		--output $OUTPUT 
	```

* **Parallelization Options**

	When running any workflow you can add the following command line options to make use of existing computing resources:
	* --local-jobs <1> : Run multiple tasks locally in parallel. Provide the max number of tasks to run at once. The default is one task running at a time.
	* --grid-jobs <0> : Run multiple tasks on a grid in parallel. Provide the max number of grid jobs to run at once. The default is zero tasks are submitted to a grid resulting in all tasks running locally.
	* --grid \<slurm> : Set the grid available on your machine. This will default to the grid found on the machine with options of slurm and sge.
	* --grid-partition \<serial_requeue> : Jobs will be submitted to the partition selected. The default partition selected is based on the default grid.

	For additional workflow options, see the [AnADAMA2](https://github.com/biobakery/anadama2) user manual.


### Standard Workflows
#### MetaWIBELE-characterize workflow
* ##### Input files for for characterization
	* clustering information for non-redudant gene catalogs using extended-fasta format, e.g. [demo_genecatalogs.clstr](https://github.com/biobakery/metawibele/examples/input/demo_genecatalogs.clstr)
	* protein seqeuences for non-redudant gene catalogs, e.g. [demo_genecatalogs.centroid.faa](https://github.com/biobakery/metawibele/examples/input/demo_genecatalogs.centroid.faa)
	* reads counts table for non-redudant gene catalogs, e.g. [demo\_genecatalogs_counts.all.tsv](https://github.com/biobakery/metawibele/examples/input/demo_genecatalogs_counts.all.tsv)
	* metadata file, e.g. [demo\_mgx_metadata.tsv](https://github.com/biobakery/metawibele/examples/input/demo_mgx_metadata.tsv)
	* sample list, e.g. [demo\_mgx_samples.tsv](https://github.com/biobakery/metawibele/examples/input/demo_mgx_samples.tsv)
	* all the above information can be specified in the `metawibele.cfg` file.

* ##### MetaWIBELE-characterize workflow
	`$ metawibele_workflow characterize --input $INPUT --output $OUTPUT`
	* Make sure the global configuration file `metawibele.cfg` is in your working directory.
	* In the command replace $INPUT with the path to the folder containing your fastq input files and $OUTPUT with the path to the folder to write output files. See the section on parallelization options to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings to run all modules. These settings will work for most data sets. However, if you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting. You can specify which modules you want to run in your own configuration file.
	* For example, `--characterization-config $myconfig_file` will modify the default settings when running the characterization modules.
	
* ##### Demo run of MetaWIBELE-characterize

	`$ metawibele_workflow characterize --input examples/input/ --output $OUTPUT_DIR`

* ##### Output files of MetaWIBELE-characterize
	**1. Annotation file**
	
	```
	familyID    annotation  feature category    method  AID
	Cluster_1   demo    study   project Shotgun NA
	Cluster_1   Cluster_1   protein_family  Denovo_clustering   CD-hit  Cluster_1__Denovo_clustering
	Cluster_1   UniRef90_A0A3E2UKI3 strong_homology UniRef90_homology   UniRef90    Cluster_1__UniRef90_homology
	Cluster_1   Faecalibacterium prausnitzii    Species Taxonomy_characterization   Taxonomy_annotation	Cluster_1__Taxonomy_characterization
	Cluster_1   UniRef90_uncharacterized    UniRef90_uncharacterized    UniRef90_characterization   UniRef90    Cluster_1__UniRef90_uncharacterized
	Cluster_1   184.02247661692778  DNA_abundance   Denovo_characterization DNA Cluster_1__DNA_abundance
	Cluster_1   0.9410658307210031  DNA_prevalence  Denovo_characterization DNA Cluster_1__DNA_prevalence
	Cluster_1   PF00408:PF02878;PF00408:PF02879;PF00408:PF02880;PF02878:PF02879;PF02878:PF02880;PF02879:PF02880 DOMINE_interaction  Denovo_characterization DOMINE  Cluster_1__DOMINE_interaction
	...
	```
	
	* File name: $OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies\_annotation.tsv
	* This is the main characterization results.
	* This file details the annotation of each protein family in the community. Protein families are groups of evolutionarily-related protein-coding sequences that often perform similar functions.
	* MetaWIBELE annotate protein family by combining global-homology similarity, local-homology similarity and non-homology based methods.
	* The annotations for each protein family coming from multiple information soruces, e.g. biochemical annotation, taxonomical annotation, ecological properties and association with host phenotypes, etc.
			
	**2. Attribute file**
	
	```
	TID AID key value
	1   Cluster_1__Denovo_clustering    repID   HSMA33LJ_27125
	2   Cluster_1__Denovo_clustering    rep_length  504
	3   Cluster_1__Denovo_clustering    cluster_size    139
	4   Cluster_1__UniRef90_homology    UniProtKB   A0A3E2UKI3_9FIRM
	5   Cluster_1__UniRef90_homology    description Phosphoglucosamine mutase
	6   Cluster_1__UniRef90_homology    organism    Faecalibacterium prausnitzii
	7   Cluster_1__UniRef90_homology    query_cov_type  high_confidence
	8   Cluster_1__UniRef90_homology    mutual_cov_type high_confidence
	9   Cluster_1__UniRef90_homology    identity    95.2
	10  Cluster_1__UniRef90_homology    query_coverage  0.9146825396825397
	11  Cluster_1__UniRef90_homology    mutual_coverage 0.9146825396825397
	12  Cluster_1__UniRef90_homology    taxa_id 1239
	13  Cluster_1__UniRef90_homology    taxa_name   Firmicutes
	...
	```
	
	* File name: $OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies\_annotation.attribute.tsv
	* This is the supplementory results for characterization.
	* Each of item is the supplemental information about the corresponding results. `AID` is the key to connect `$BASENAME_proteinfamilies_annotation.tsv` with `$BASENAME_proteinfamilies_annotation.attribute.tsv`.
	
	**3. Taxonomic file**
	
	```
	familyID    study   map_type    query_type  mutual_type identity    query_coverage  mutual_coverage detail  Tax TaxID   Rep_Tax Rep_TaxID   organism    UniProtKB   unirefID    note    msp_name    msp_taxa_name   msp_taxa_id MSP_Tax MSP_TaxID   MSP_Rep_Tax MSP_Rep_TaxID   taxa_id taxa_name   taxa_rank   taxa_lineage
	Cluster_1   demo    UniRef90_uncharacterized    high_confidence high_confidence 95.2    0.9146825396825397  0.9146825396825397  Phosphoglucosamine mutase   Firmicutes  1239    Faecalibacterium prausnitzii    853 Faecalibacterium prausnitzii    A0A3E2UKI3_9FIRM    UniRef90_A0A3E2UKI3 good    msp_unknown NA  NA  Faecalibacterium prausnitzii    853 Faecalibacterium prausnitzii    853 853 Faecalibacterium prausnitzii    Species k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii
	Cluster_10  demo    UniRef90_uncharacterized    high_confidence high_confidence 93.3    0.9782608695652174  0.9782608695652174  MafF    Bacteria    2   Escherichia coli    562 Escherichia coli    B8QUG6_ECOLX    UniRef90_B8QUG6 good    msp_unknown NA  NA  NA  NA  Escherichia coli    562 562 Escherichia coli    Species k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli
	Cluster_100 demo    UniRef90_characterized  high_confidence high_confidence 95.8    1.0 1.0 Uncharacterized protein Clostridiales   186802  Faecalibacterium prausnitzii M21/2  411485  Faecalibacterium prausnitzii M21/2  A8SEK6_9FIRM    UniRef90_A8SEK6 good    msp_unknown NA  NA  Faecalibacterium prausnitzii M21/2  411485  Faecalibacterium prausnitzii M21/2  411485  411485  Faecalibacterium prausnitzii M21/2  Terminal    k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|t__Faecalibacterium_prausnitzii_M21/2
	Cluster_1000    demo    UniRef90_characterized  high_confidence high_confidence 99.8    1.0 1.0 FeS assembly protein SufD   Bacteroidaceae  815 Bacteroides fragilis HMW 615    1073387 Bacteroides fragilis HMW 615    K1GB77_BACFG    UniRef90_K1GB77 good    msp_02  Bacteroides fragilis    817 Bacteroides 816 Bacteroides fragilis HMW 615    1073387 1073387 Bacteroides fragilis HMW 615    Terminal    k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_fragilis|t__Bacteroides_fragilis_HMW_615
	Cluster_10000   demo    UniRef90_characterized  high_confidence high_confidence 90.9    0.900355871886121   0.900355871886121   Methanol dehydrogenase  Faecalibacterium    216851  Faecalibacterium prausnitzii    853 Faecalibacterium prausnitzii    A0A329U1M8_9FIRM    UniRef90_A0A329U1M8 good    msp_unknown NA  NA  NA  NA  Faecalibacterium prausnitzii    853 853 Faecalibacterium prausnitzii    Species k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii
	...
	```
	
	* File name: $OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies\_annotation.MSPminer\_taxonomy.tsv
	* These files report the detailed information about taxonomic annotation.
	* `$BASENAME_proteinfamilies_annotation.MSPminer_taxonomy.all.tsv` shows the taxonomic information at each of taxonomic level.
	
	
	**4. Abundance file**
	
	```
	ID      CSM5FZ3N_P      CSM5FZ3R_P      CSM5FZ3T_P      CSM5FZ3V_P      CSM5FZ3X_P      CSM5FZ3Z_P      CSM5FZ42_P      CSM5F
	Cluster_1       399.304 15.8845 0       171.157 1.41663 0.275544        1.16031 0       0       8.46607 321.148 60.9853 347.5
	Cluster_10      11.1512 8.81212 9.00341 5.77603 0       3.019   4.62291 5.2878  2.12412 7.00065 8.10752 27.1621 47.862  39.04
	Cluster_100     54.6202 26.7423 0       13.5309 8.26367 0       2.21515 11.261  0       44.7264 12.0862 65.9435 284.451 165.2
	Cluster_1000    21.6101 11.235  134.379 157.494 244.854 361.196 93.9497 126.383 269.297 0       280.881 13.575  11.6934 20.65
	Cluster_10000   62.6742 4.6883  0       48.853  0       0       0       0       0       0       72.9965 2.8902  20.9542 16.77
	Cluster_10001   2.42532 7.18719 0       0.628126        0       0       0       0       0       12.5615 7.93502 2.21535 8.351
	...
	```
	
	* File name: $OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies\_relab.tsv
	* This is the relative abundance per protein family across samples.
	* Protein family abundance is reported in copies per million (CPM) units, which is "total sum scaling (TSS)"-style normalization: each sample is constrained to sum to 1 million. First, each protein family is normalized to RPK (reads per kilobase) units for gene length normalization; RPK units reflect relative gene (or transcript) copy number in the community. Then, RPK values are further sum-normalized (CPM) to adjust for differences in sequencing depth across samples. Further information can refer to the normalization approach in [HUMAnN](https://github.com/biobakery/humann). 
	
	**5. Clusting information for protein families**
	
	```
	>HSMA33LJ_27125;Cluster_1;length=504;size=139;cluster=1
	HSMA33LJ_27125
	HSM5MD82_P_40510
	HSM6XRQ8_06389
	ESM5GEYY_P_120902
	HSM5MD5F_P_64814
	MSM9VZIU_90100
	...
	```
	
	* File name: $OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies.clstr
	* This is the clustering information for protein families, formatted using extention-fasta style based on the version of CD-hit clustering file.
	
	**6. Sequences of protein families**  
	
	* File name: $OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies.centroid.faa
	* This file is the protein sequences for representatives of protein families.

	**7. Intermediate output files**
	
	* Clustering results
		* MetaWIBELE clusters all representative sequences of gene catalogs into protein families. 
		* All intermediate results are in the folder `$OUTPUT_DIR/characterization/clustering`.
		
	* Protein-family search results
		* MetaWIBELE queries each sequence in protein families against UniRef90 database by performing protein-level search.
		* All intermediate results are in the folder `$OUTPUT_DIR/characterization/protein_family_annotation`.
	
	* Domain-motif annotation results
		* MetaWIBELE uses a local-homology approach (domain/motif search) to characterize the secondary structures of protein families.
		* All intermediate results are in the folder `$OUTPUT_DIR/characterization/domain_motif_annotation`.
	
	* Abundance-based annotation results
		* MetaWIBELE implements a non-homology based strategy compromising (i) taxonomic annotation with phylogenetic binning, (ii) abundance profiling for protein families and (iii) association with host phenotypes based on differential abundance. 
		* All intermediate results are in the folder `$OUTPUT_DIR/characterization/abundance_annotation`.
	

#### MetaWIBELE-prioritize workflow
* ##### Input files for prioritiztion
	* anntation file produced by MetaWIBELE-characterize workflow:$OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies\_annotation.tsv
	* attribute file produced by MetaWIBELE-characterize workflow: $OUTPUT_DIR/characterization/$BASENAME\_proteinfamilies\_annotation.attribute.tsv

* ##### MetaWIBELE-prioritize workflow

	`$ metawibele_workflow prioritize --input $INPUT --output $OUTPUT`
	* In the command replace $INPUT with the path to the folder containing your fastq input files and $OUTPUT with the path to the folder to write output files. See the section on parallelization options to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings for all main tool subtasks. These settings will work for most data sets. However, if you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting. Then apply these settings by using options for each task. You can specify your own configuration file.
	* For example, `--prioritization-config $myconfig_file` will modify the default settings when running the prioritization tasks.
	
* ##### Demo run of MetaWIBELE-prioritize

	`$ metawibele_workflow prioritize --input examples/characterization/ --output $OUTPUT\_DIR/`

* ##### Output files of MetaWIBELE-prioritize
	**1. unsupervised prioritization**
	
	```
	familyID    DNA_nonIBD_abundance__value DNA_nonIBD_abundance__percentile    DNA_nonIBD_prevalence__value    DNA_nonIBD_prevalence__percentile   priority_score
	Cluster_24570   2694.0590678779345  0.9999793230362054  1.0 1.0 0.9999896614112173
	Cluster_41147   2431.225870892018   0.9999586460724107  1.0 1.0 0.9999793226086595
	Cluster_22422   1336.999313239437   0.9998966151810268  1.0 1.0 0.9999483049182701
	Cluster_40049   1109.6718042394361  0.9998759382172322  1.0 1.0 0.9999379652605459
	Cluster_29449   803.4913150469486   0.9997311994706697  0.9976525821596244  0.9993898213934824  0.999560481284518
	Cluster_21419   383.6281626291082   0.9983458428964291  1.0 1.0 0.9991722368230451
	...
	```
	
	* File name: 
	$OUTPUT\_DIR/prioritization/$BASENAME\_unsupervised\_prioritization.rank.tsv
	* These are the results of unsupervised prioritization based on ecological properties. Each of protein family has a numeric priority score.
	* `$BASENAME_unsupervised_prioritization.rank.tsv` is the overall ranking for all protein families.

	
	**2. supervised prioritization: numeric ranking**
	
	```
	familyID    DNA_within_phenotype_abundance__value   DNA_within_phenotype_abundance__percentile  DNA_within_phenotype_prevalence__value  DNA_within_phenotype_prevalence__percentile MaAsLin2_DA__coef__value    MaAsLin2_DA__coef__percentile   MaAsLin2_DA__qvalue__value  MaAsLin2_DA__qvalue__percentile priority_score
	Cluster_14393|CD.dysbiosis_vs_CD.non_dysbiosis  844.0252184037556   0.9995526838966203  0.971830985915493   0.986622572543949   -480.492828810768   0.9995526838966203  4.31866766640033e-11    1.0 0.9963995495816218
	Cluster_47254|CD.dysbiosis_vs_CD.non_dysbiosis  718.2714984741792   0.9992047713717693  0.9741784037558685  0.9884128602332347  -392.524817430569   0.9993538767395627  4.87065570823145e-09    0.9937619603847205  0.9951628678096055
	Cluster_53|CD.dysbiosis_vs_CD.non_dysbiosis 357.3594103568074   0.9953280318091451  0.9954233409610984  0.9983340378446925  -143.260405479191   0.9928429423459244  4.93693768852676e-09    0.9936625493948356  0.9950374594865485
	Cluster_7|CD.dysbiosis_vs_CD.non_dysbiosis  438.19562671167034  0.9975646123260438  0.9741784037558685  0.9884128602332347  -250.126179947757   0.9979622266401591  3.60646971240824e-09    0.9947063647886274  0.9946467984814091
	Cluster_14|CD.dysbiosis_vs_CD.non_dysbiosis 428.82273066590415  0.9973658051689861  0.9624413145539906  0.9785911430489594  -256.077370752594   0.9982107355864811  3.00504043491864e-09    0.9950294505057534  0.9922342227131524
	Cluster_42|CD.dysbiosis_vs_CD.non_dysbiosis 437.4689617025174   0.9975149105367793  0.9765258215962441  0.9895317900390382  -220.208838202021   0.9974652087475149  1.33219876365297e-07    0.979198250366578   0.9908703385076423
	...
	```
	
	* File name: 
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.tsv
	* These are the results of supervised prioritization by combing ecological properties and assoiciation with host phenotypes. Each of protein family has a numeric priority score.
	* `$BASENAME_supervised_prioritization.rank.tsv` is the overall ranking for all protein families.


	**3. supervised prioritization: binary filtering**
	
	***3.1 Select interested subset annotated with at least one of specific biochemical annotations***
	 
	Setting `prioritization.cfg` as following:
	
	```
	##Binary filtering for selection subset
	[filtering]
	# Filter for interested functional vignettes type [Choices: pilin | secreted_system | other user defined | none]
	vignettes = none

	# significant associations for filtering: [required] required item, [optional] optional item, [none] ignoring
	MaAsLin2_DA-sig = none

	# biochemical annotation for filtering: [required] required item, [optional] optional item, [none] ignoring
	ExpAtlas_interaction = required
	DOMINE_interaction = none
	SIFTS_interaction = none
	Denovo_signaling = optional
	Denovo_transmembrane = optional
	PSORTb_extracellular = optional
	PSORTb_cellWall = optional
	PSORTb_outerMembrane = optional
	UniRef90_extracellular = optional
	UniRef90_signaling = optional
	UniRef90_transmembrane = optional
	UniRef90_cellWall = optional
	UniRef90_outerMembrane = optional
	UniRef90_PfamDomain = none
	InterProScan_PfamDomain = none
	InterProScan_SUPERFAMILY = none
	InterProScan_ProSiteProfiles = none 
	InterProScan_ProSitePatterns = none
	InterProScan_Gene3D = none
	InterProScan_PANTHER = none
	InterProScan_TIGRFAM = none
	InterProScan_SFLD = none
	InterProScan_ProDom = none
	InterProScan_Hamap = none
	InterProScan_SMART = none
	InterProScan_CDD = none
	InterProScan_PRINTS = none
	InterProScan_PIRSF = none
	InterProScan_MobiDBLite = none
	```
	
	* Re-run prioritization workflow for filtering:
		* `$ metawibele_workflow prioritize --prioritization-config prioritization.cfg --bypass-mandatory --selected-output demo_prioritized.selected.tsv --output $OUTPUT\_DIR/`
	* Output file name: $OUTPUT_DIR/prioritization/demo_prioritized.selected.tsv
	* This file is the results of supervised filtering of protein families based on biochemical annotations.
	* These settings require that each of prioritized protein family should 1) be annotated to domain-domain interaction with host, and 2) have at least one of the following features: signaling, extracellular, cellWall, outerMembrane, transmembrane 
	
    ***3.2. Select interested subset annotated with multiple specific biochemical annotations simultaneously***
	 
	 Setting `prioritization.cfg` as following:
	
	```
	##Binary filtering for selection subset
	[filtering]
	# Filter for interested functional vignettes type [Choices: pilin | secreted_system | other user defined | none]
	vignettes = none

	# significant associations for filtering: [required] required item, [optional] optional item, [none] ignoring
	MaAsLin2_DA-sig = required

	# biochemical annotation for filtering: [required] required item, [optional] optional item, [none] ignoring
	ExpAtlas_interaction = required
	DOMINE_interaction = none
	SIFTS_interaction = none
	Denovo_signaling = required
	Denovo_transmembrane = optional
	PSORTb_extracellular = optional
	PSORTb_cellWall = optional
	PSORTb_outerMembrane = optional
	UniRef90_extracellular = optional
	UniRef90_signaling = none
	UniRef90_transmembrane = optional
	UniRef90_cellWall = optional
	UniRef90_outerMembrane = optional
	UniRef90_PfamDomain = none
	InterProScan_PfamDomain = none
	InterProScan_SUPERFAMILY = none
	InterProScan_ProSiteProfiles = none 
	InterProScan_ProSitePatterns = none
	InterProScan_Gene3D = none
	InterProScan_PANTHER = none
	InterProScan_TIGRFAM = none
	InterProScan_SFLD = none
	InterProScan_ProDom = none
	InterProScan_Hamap = none
	InterProScan_SMART = none
	InterProScan_CDD = none
	InterProScan_PRINTS = none
	InterProScan_PIRSF = none
	InterProScan_MobiDBLite = none
	```
	
	 * Re-run prioritization workflow for filtering:
		* `$ metawibele_workflow prioritize --prioritization-config prioritization.cfg --bypass-mandatory --selected-output demo_prioritized.selected.tsv --output $OUTPUT\_DIR/`
	 * Output file name: $OUTPUT_DIR/prioritization demo_prioritized.selected.tsv
	 * This file is the results of supervised filtering of protein families based on biochemical annotations.
	 * These settings require that each of prioritized protein family should 1) significantly associated with the main phenotype, 2) be annotated to domain-domain interaction with host, 3) predicted as signal peptides, and 4) have at least one of the following features: extracellular, cellWall, outerMembrane, transmembrane 
	
	***3.3 Select interested subset based on specific functions***
	
	Setting `prioritization.cfg` as following:
	
	```
	[filtering]
	# Filter for interested functional vignettes type [Choices: pilin | secreted_system | other user defined | none]
	vignettes = pilin

	# Filter for significant associations: [required] required item, [optional] optional item, [none] ignoring [ Default: required ]
	MaAsLin2_DA-sig = none

	# biochemical annotation for filtering: [required] required item, [optional] optional item, [none] ignoring
	ExpAtlas_interaction = optional
	DOMINE_interaction = optional
	SIFTS_interaction = optional
	Denovo_signaling = optional
	Denovo_transmembrane = optional
	PSORTb_extracellular = optional
	PSORTb_cellWall = optional
	PSORTb_outerMembrane = optional
	UniRef90_extracellular = optional
	UniRef90_signaling = optional
	UniRef90_transmembrane = optional
	UniRef90_cellWall = optional
	UniRef90_outerMembrane = optional
	UniRef90_PfamDomain = optional
	InterProScan_PfamDomain = optional
	InterProScan_SUPERFAMILY = optional
	InterProScan_ProSiteProfiles = optional 
	InterProScan_ProSitePatterns = optional
	InterProScan_Gene3D = optional
	InterProScan_PANTHER = optional
	InterProScan_TIGRFAM = optional
	InterProScan_SFLD = optional
	InterProScan_ProDom = optional
	InterProScan_Hamap = optional
	InterProScan_SMART = optional
	InterProScan_CDD = optional
	InterProScan_PRINTS = optional
	InterProScan_PIRSF = optional
	InterProScan_MobiDBLite = optional
	```
	
	* Re-run prioritization workflow for filtering:
		* `$ metawibele_workflow prioritize --prioritization-config prioritization.cfg --bypass-mandatory --vignette-config my_vignette_function_file  --selected-output demo_prioritized_pilin.tsv --output $OUTPUT\_DIR/`
	* Provide your own vignette function file for filtering specific functions.
	* Output file name: $OUTPUT_DIR/prioritization/demo_prioritized_pilin.tsv
	* This file is the results of supervised filtering of protein families based on pilin related functions.
	
	
	**4. finalized prioritization**
	
	```
	TID familyID evidence    value   rank    description note
	1   Cluster_14393   DNA_within_phenotype_abundance  844.0252184037556   0.9995526838966203	ranking based on single evidence  CD.dysbiosis_vs_CD.non_dysbiosis
	2   Cluster_14393   DNA_within_phenotype_prevalence 0.971830985915493   0.986622572543949	ranking based on single evidence   CD.dysbiosis_vs_CD.non_dysbiosis
	3   Cluster_14393   MaAsLin2_DA__coef   -480.492828810768   0.9995526838966203	ranking based on single evidence  CD.dysbiosis_vs_CD.non_dysbiosis
	4   Cluster_14393   MaAsLin2_DA__qvalue 4.31866766640033e-11    1.0	ranking based on single evidence CD.dysbiosis_vs_CD.non_dysbiosis
	5   Cluster_14393   priority_score  0.9963995495816218  0.9963995495816218  meta ranking based on multiple evidences    CD.dysbiosis_vs_CD.non_dysbiosis
	...
	```
	
	* File name: 
		$OUTPUT\_DIR/prioritization/$BASENAME\_unsupervised\_prioritization.rank.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.selected.table.tsv
	* These file are formated the prioritization results in the same way.

***


## Guides to MetaWIBELE Utilities

### Preprocessing sequencing reads to build gene catalogs
A utility workflow in MetaWIBELE package for preprocessing metagenomes reads, used for (i) metagenomic assembly, (ii) open reading frame prediction, (iii) non-redundant gene catalogs construction and (iv) gene abundance estimation.

#### Preprocessing workflow
```
usage: preprocess.py [-h] [--version] [--threads THREADS]
                     [--extension-paired EXTENSION_PAIRED]
                     [--sample-list SAMPLE_LIST]
                     [--gene-call-type {prokka,prodigal,both}]
                     [--extension {.fastq.gz,.fastq}] [--bypass-assembly]
                     [--bypass-gene-calling] [--bypass-gene-catalog]
                     [--output-basename OUTPUT_BASENAME] -o OUTPUT [-i INPUT]
                     [--config CONFIG] [--local-jobs JOBS]
                     [--grid-jobs GRID_JOBS] [--grid GRID]
                     [--grid-partition GRID_PARTITION]
                     [--grid-benchmark {on,off}] [--grid-options GRID_OPTIONS]
                     [--grid-environment GRID_ENVIRONMENT]
                     [--grid-scratch GRID_SCRATCH] [--dry-run]
                     [--skip-nothing] [--quit-early] [--until-task UNTIL_TASK]
                     [--exclude-task EXCLUDE_TASK] [--target TARGET]
                     [--exclude-target EXCLUDE_TARGET]
                     [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

A workflow to preprocess shotgun sequencing reads of metagenomes with tasks of metagenomic assembly, gene calling, building gene catalogs and generating gene abundance for each sample.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --threads THREADS     number of threads/cores for each task to use
  --extension-paired EXTENSION_PAIRED
                        provide the extension for paired fastq files using comma to seperate, e.g. .R1.fastq.gz,.R2.fastq.gz | .R1.fastq,.R2.fastq
  --sample-list SAMPLE_LIST
                        sample list file
  --gene-call-type {prokka,prodigal,both}
                        specify which type of gene calls will be used
                        [default: both]
  --extension {.fastq.gz,.fastq}
                        provide the extension for all fastq files
                        [default: .fastq.gz]
  --bypass-assembly     do not run assembly
  --bypass-gene-calling
                        do not call ORFs
  --bypass-gene-catalog
                        do not build gene catalogs
  --output-basename OUTPUT_BASENAME
                        provide the basename for output files
  -o OUTPUT, --output OUTPUT
                        Write output to this directory
  -i INPUT, --input INPUT
                        Find inputs in this directory 
                        [default: /srv/export/hutlab11/share_root/users/yancong/demo]
  --config CONFIG       Find workflow configuration in this folder 
                        [default: only use command line options]
  --local-jobs JOBS     Number of tasks to execute in parallel locally 
                        [default: 1]
  --grid-jobs GRID_JOBS
                        Number of tasks to execute in parallel on the grid 
                        [default: 0]
  --grid GRID           Run gridable tasks on this grid type 
                        [default: slurm]
  --grid-partition GRID_PARTITION
                        Partition/queue used for gridable tasks.
                        Provide a single partition or a comma-delimited list
                        of short/long partitions with a cutoff.
                        [default: serial_requeue,shared,240]
  --grid-benchmark {on,off}
                        Benchmark gridable tasks 
                        [default: on]
  --grid-options GRID_OPTIONS
                        Grid specific options that will be applied to each grid task
  --grid-environment GRID_ENVIRONMENT
                        Commands that will be run before each grid task to set up environment
  --grid-scratch GRID_SCRATCH
                        The folder to write intermediate scratch files for grid jobs
  --dry-run             Print tasks to be run but don't execute their actions 
  --skip-nothing        Run all tasks. Rerun tasks that have already been run.
  --quit-early          Stop if a task fails. By default,
                        all tasks (except sub-tasks of failed tasks) will run.
  --until-task UNTIL_TASK
                        Stop after running this task. Use task name or number.
  --exclude-task EXCLUDE_TASK
                        Don't run these tasks. Add multiple times to append.
  --target TARGET       Only run tasks that generate these targets.
                        Add multiple times to append.
                        Patterns with ? and * are allowed.
  --exclude-target EXCLUDE_TARGET
                        Don't run tasks that generate these targets.
                        Add multiple times to append.
                        Patterns with ? and * are allowed.
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the level of output for the log 
                        [default: INFO]
```

* `--input`: the input directory where a set of fastq (or fastq.gz) files (single-end or paired-end) passing through QC are stored. The files are expected to be named `$SAMPLE.paired_R1.gz`, `$SAMPLE.paired_R2.gz`, `$SAMPLE.orphan_R1.gz` and `$SAMPLE.orphan_R2.gz` where `$SAMPLE` is the sample name or identifier corresponding to the sequences. `$SAMPLE` can contain any characters except spaces or periods.
* `--extension-paired` indicates the extension for paired fastq files using comma to seperate. It should be specified as ".R1.fastq.gz,.R2.fastq.gz" if the paired fastq files are `$SAMPLE.R1.fastq.gz` and `$SAMPLE.R2.fastq.gz`  
* `--extension` indicates the extension for all fastq files. It should be specified as ".fastq.gz" if the fastq files are `$SAMPLE.fastq.gz` 
* `--output`: the ouput directory. 

#### Input files of preprocessing workflow
* Make sure the `metawibele.cfg` in your working directory.
* QC'ed shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format)
* See the section on parallelization options to optimize the workflow run based on your computing resources. 
* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting, you can change the parameter settings.
* For example, `--extension-paried "$R1_suffix,$R2_suffix"`, `--extension "$fastq_suffix"` (what are the follwong part after `$SAMPLE` in the input file names) will modify the default settings when running the assembly task.

#### Demo run of preprocessing workflow

`$ preprocessing_workflow --input examples/raw_reads/ --output preprocessing/ --extension-paired "_R1.fastq.gz,_R2.fastq.gz" --extension ".fastq.gz" --sample-list examples/raw_reads/sample.txt --local-jobs 5 --output $OUTPUT_DIR/`

#### Output files of preprocessing workflow
**1. assembly results**
	
* `$OUTPUT_DIR/assembly/$BASENMAE_contig_sequence.fasta`: contig sequences
* The assembly outputs for each of sample are in `$OUTPUT_DIR/assembly/` folder.
	
**2. gene-calling results**
	
* `$OUTPUT_DIR/$BASENMAE_gene_info.tsv`: all gene calls information
* `$OUTPUT_DIR/$BASENMAE_combined_gene_protein_coding.sorted.fna`: nucleotide sequences for all ORFs sorted.
* `$OUTPUT_DIR/$BASENMAE_combined_protein.sorted.faa`: protein sequences for all ORFs sorted by gene length.
* `$OUTPUT_DIR/$BASENMAE_combined_gene_protein_coding.complete.sorted.fna`: nucleotide sequences for all complete ORFs sorted by gene length.
* `$OUTPUT_DIR/$BASENMAE_combined_protein.complete.sorted.faa`: protein sequences for all complete ORFs sorted by gene length.
* The gene-calling outputs from prodigal are in `$OUTPUT_DIR/gene_calls` folder. 
* The gene-annotation outputs from prokka are in `$OUTPUT_DIR/gene_annotation` folder.
	
**3. gene catalogs**
	
* `$OUTPUT_DIR/$BASENMAE_genecatalogs.clstr`: clustering information for non-redudant gene catalogs
* `$OUTPUT_DIR/$BASENMAE_genecatalogs.centroid.fna`: nucleotide sequences of representatives for gene catalogs.
* `$OUTPUT_DIR/$BASENMAE_genecatalogs.centroid.faa`: protein sequences of representatives for gene catalogs.
* `$OUTPUT_DIR/$BASENMAE_genecatalogs_counts.all.tsv`: reads counts of gene catalogs across samples.
* All mapping outputs for each sample are in `$OUTPUT_DIR/mapping` folder. 
	
***


## Download MetaWIBELE resources
### Information of gene catalogs
* [HMP2\_contig_sequence.fasta.tar.gz]() (7.5 GB): contig sequences
* [HMP2\_gene_info.tsv.tar.gz]() (2.3 GB): information of gene calling
* [HMP2_genecatalogs.clstr.tar.gz]() (279 MB): clustering information for gene catalogs
* [HMP2_genecatalogs.centroid.fna.tar.gz]() (554 MB): nucleotide sequences of centroids for gene catalogs
* [HMP2_genecatalogs.centroid.faa.tar.gz]() (335 MB): protein sequences of centroids for gene catalogs
* [HMP2\_genecatalogs_CPM.tsv.tar.gz]() (2.0 GB): DNA relative abundance of gene catalogs

### Characterization of protein families 
* [HMP2_proteinfamilies.clstr.tar.gz]() (29 MB): clustering information for protein families
* [HMP2\_proteinfamilies.centroid.faa.tar.gz]() (270 MB): protein sequences of centroids for protein families
* [HMP2\_proteinfamilies_CPM.tsv.tar.gz]() (1.5 GB): DNA relative abundance of protein families
* [HMP2\_proteinfamilies_annotation.tsv.tar.gz]() (1.8 GB): main annotations of protein families
* [HMP2\_proteinfamilies_annotation.attribute.tsv.tar.gz]() (6.4 GB): attributes of annotation types

### Prioritization of protein families
* [HMP2\_unsupervised_prioritization.rank.table.tsv.tar.gz]() (14 MB): prioritization based on overall abundance and prevalence
* [HMP2\_supervised_prioritization.rank.table.tsv.tar.gz]() (13 MB): prioritization based on assocation with phenotypes and abundance information 
# MetaWIBELE User Manual

**MetaWIBELE** (**W**orkflow to **I**dentify novel **B**ioactive **Ele**ments in the microbiome) is a workflow to efficiently and systematically identify and prioritize potentially bioactive (and often uncharacterized) gene products in microbial communities. It prioritizes candidate gene products from assembled metagenomes using a combination of sequence homology, secondary-structure-based functional annotations, phylogenetic binning, ecological distribution, and association with environmental parameters or phenotypes to target candidate bioactives.


## Citing MetaWIBELE

**A manuscript describing MetaWIBELE is currently in prep:**

Identifying Novel Bioactive Microbial Gene Products in Inflammatory Bowel Disease

**And feel free to link to MetaWIBELE in your Methods:**

[http://huttenhower.sph.harvard.edu/metawibele](http://huttenhower.sph.harvard.edu/metawibele)


**For additional information, read the [MetaWIBELE Tutorial](https://github.com/biobakery/biobakery/wiki/metawibele).**


If you have questions, please direct it to [the MetaWIBELE channel](https://forum.biobakery.org/c/Microbial-community-profiling/metawibele) of the bioBakery Support Forum.

***

## Contents ##

* [Workflow](#workflow)
	* [Workflow by bypass mode](#workflow-by-bypass-mode) 
* [Install MetaWIBELE](#install-metawibele)
    * [Requirements](#requirements)
    * [Installation](#installation)
    	* [Download MetaWIBELE](#download-metawibele)
    	* [Install MetaWIBELE](#install-metawibele)
    	* [Install databases](#install-databases)
    		* [UniRef database](#uniref-database)
    		* [Domain database](#domain-database)
    	* [Prepare configuration files](#prepare-configuration-files)
			* [Prepare global configuration file](#prepare-global-configuration-file)
    		* [Prepare local configuration file](#prepare-local-configuration-file)
    		* [Prepare vignette configuration file](#prepare-vignette-configuration-file)
* [Quick-start Guide](#quick-start-guide)
    * [How to run](#how-to-run)
    * [Standard Workflows](#standard-workflows)
    	* [MetaWIBELE-characterize](#metawibele-characterize)
    		* [Input files for characterization](#input-files-for-characterization)
    		* [MetaWIBELE-characterize workflow](#metawibele-characterize-workflow)
    		* [Demo run of MetaWIBELE-characterize](#demo-run-of-metawibele-characterize)
    		* [Output files of MetaWIBELE-characterize](#output-files-of-metawibele-characterize)
    	* [MetaWIBELE-prioritize](#metawibele-prioritize)
    		* [Input files for prioritization](#input-files-for-prioritization)
    		* [MetaWIBELE-prioritize workflow](#metawibele-prioritize-workflow)
    		* [Demo run of MetaWIBELE-prioritize](#demo-run-of-metawibele-prioritize)
    		* [Output files of MetaWIBELE-prioritize](#output-files-of-metawibele-prioritize)
* [Guides to MetaWIBELE Utilities](#guides-to-metawibele-utilities)
	* [Preprocessing sequencing reads to build gene families](#preprocessing-sequencing-reads-to-build-gene-families)
		* [Preprocessing workflow](#preprocessing-workflow)
		* [Input files for preprocessing workflow](#input-files-preprocessing-workflow)
		* [Demo run of preprocessing workflow](#demo-run-of-preprocessing-workflow) 
		* [Output files of preprocessing workflow](#output-files-of-preprocessing-workflow)
    	
***


## Workflow
### Description
MetaWIBELE identifies, prioritizes, and putatively annotates potentially bioactive gene products in host- and non-host-associated microbial communities and human disease phenotypes. This software starts with gene products from assembled metagenomes and combines sequence homology, secondary-structure-based functional annotations, taxonomic profiling, ecological distribution, and environmental or phenotypic association to identify candidate bioactives.


### Workflow by bypass mode
There are multiple bypass options that will allow you to adjust the standard workflow.

Bypass options:

* --bypass-global-homology
	* do not annotate protein families based on global homology information
* --bypass-domain-motif
	* do not annotate protein families based on domain/motif information
* --bypass-interproscan
	* do not annotate protein families based on Interproscan for domain/motif characterization
* --bypass-pfam\_to\_go
	* do not annotate protein families based on Pfam2GO for domain/motif characterization
* --bypass-domine
	* do not annotate protein families based on DOMINE database for domain/motif characterization
* --bypass-sifts
	* do not annotate protein families based on SIFTS database for domain/motif characterization
* --bypass-expatlas
	* do not annotate protein families based on Expression Atlas database for domain/motif characterization
* --bypass-psortb
	* do not annotate protein families based on PSORTb for domain/motif characterization	
* --bypass-abundance
	* do not annotate protein families based on abundance information
* --bypass-integration
	* do not integrate annotations for protein families
* --bypass-mandatory
	* do not prioritize protein families to calculate continuous priority scores
* --bypass-optional
	* do not prioritize protein families based on selecting our for interested annotations (optional prioritization)


## Install MetaWIBELE
### Requirements

1. [Python](https://www.python.org/) (version >= 3.7, required numpy, pandas, scipy packages)
2. [AnADAMA2](https://huttenhower.sph.harvard.edu/anadama2) (version >= 0.7.4-devel)
3. [CD-hit](http://weizhongli-lab.org/cd-hit/) (version >= 4.7)
4. [Diamond](http://www.diamondsearch.org/index.php) (version >= 0.9.24)
5. [MSPminer](https://www.enterome.com/downloads/) (version >= 2; licensed software)
6. [MaAsLin2](https://huttenhower.sph.harvard.edu/maaslin2) (version >= 1.1.2) (only required if using MaAsLin2 to associate with environmental parameters or phenotypes)
7. [Interproscan](https://github.com/ebi-pf-team/interproscan/wiki) (version >= 5.31-70) (only required if using Interproscan to annotate domains and motifs; installing with activating Phobius/SignalP/TMHMM analyses)
8. [Signalp](http://www.cbs.dtu.dk/services/SignalP-4.1/) (version >= 4.1) (only required if using Signalp to annotate signal peptides integrated in interproscan; licensed software)
9. [TMHMM](http://www.cbs.dtu.dk/services/TMHMM/) (version >= 2.0) (only required if using TMHMM to annotate transmembrane proteins integrated by interproscan; licensed software)
10. [Phobius](http://phobius.sbc.su.se/) (version >= 1.01) (only required if using Phobius to annotate both signal peptides and transmembrane proteins integrated by interproscan; licensed software)
11. [PSORTb](https://psort.org/documentation/index.html) (version >= 3.0) (only required if using PSORTb to predict subcellular localization; licensed software)
12. **Optional**: only required if using MetaWIBELE utility to preprocess metagenomic sequencing reads
	* [MEGAHIT](https://github.com/voutcn/megahit) (version >= 1.1.3) 
	* [Prokka](https://github.com/tseemann/prokka) (version >= 1.14-dev; recommend to not set '-c' parameter when running prodigal with metagenome mode)
	* [Prodigal](https://github.com/hyattpd/Prodigal) (version >= 2.6)
	* [USEARCH](http://www.drive5.com/usearch/) (version >= 9.0.2132; licensed software)
	* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.3.2)
	* [SAMtools](https://github.com/samtools/) (version >= 1.9)
	* [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) (version >= 1.6.2)

**Note:** Please install the required software in a location in your `$PATH`. If you always run with gene families (non-redundant gene catalogs), the optional softwares are not required. Also if you always run with one or more bypass options (for information on bypass options, see optional arguments to the section [Workflow by bypass mode](#workflow-by-bypass-mode)), the software required for the steps you bypass does not need to be installed.


### Installation
#### Download MetaWIBELE
You can download the latest MetaWIBELE release or the development version. The source contains example files.

Option 1: Latest Release (Recommended)

* download [metawibele-master.zip](https://github.com/biobakery/metawibele/archive/master.zip) and unpack the latest release of MetaWIBELE

Option 2: Development Version

* Create a clone of the Git repository 
	* `$ git clone https://github.com/biobakery/metawibele.git`
* You can always update to the latest version of the repository with:
	* `$ git pull --update`


#### Install MetaWIBELE
You only need to do **any one** of the following options to install the MetaWIBELE package.


Option 1: Installing with pip

* `$ pip install metawibele`
* If you do not have write permissions to `/usr/lib/`, then add the option --user to the install command. This will install the python package into subdirectories of `~/.local/`. Please note when using the --user install option on some platforms, you might need to add `~/.local/bin/` to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `metawibele: command not found` when trying to run MetaWIBELE after installing with the --user option.

Option 2: Installing with docker

* `$ docker pull biobakery/metawibele`
* This docker image includes most of the dependent software packages.
* Large software packages and those with licenses are not included in this image:
	* Softwares with the license : mspimer, signalp, TMHMM, phobius, psortb
	* Softwares with large size: interproscan
	* Users should review the license terms and install these packages manually. 
	
Option 3: Installing with conda

* `$ conda install -c biobakery metawibele`


Option 4: Installing from source

* Move to the MetaWIBELE directory
	* `$ cd $MetaWIBELE_PATH`

* Install MetaWIBELE package
	* `$ python setup.py install`
	* If you do not have write permissions to `/usr/lib/`, then add the option --user to the install command. This will install the python package into subdirectories of `~/.local/`. Please note when using the --user install option on some platforms, you might need to add `~/.local/bin/` to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `metawibele: command not found` when trying to run MetaWIBELE after installing with the --user option. 
	* Similarly, you can also specify the installation directory using --prefix option which will install the python package into the directory `$YOUR_INSTALL_DIR` that is a directory on PYTHONPATH and which Python reads ".pth" files from. You might need to add `$YOUR_INSTALL_DIR/bin` to your `$PATH` as it might not be included by default.
 

#### Install databases
To run metawibele, you need to install the dependent databases: 1) uniref databases (**required**); 2) domain databases (optional).

##### UniRef database
UniRef databases are **required** if you will use MetaWIBELE to do global-homology based annotation and taxonomic annotation. You can use **any one** of the following options to install these databases.

Option 1: Download uniref databases (Recommended)

* We have built the dependent UniRef databases based on UniProt/UniRef 2019_01 sequences and annotations. You can download and uncompress these databases (both sequences and annotations) and provide `$UNIREF_LOCATION` as the location to install the databases.
	
	* UniRef90 sequence file (20 GB): 
		* If you are using Diamond v0.9.24, just download and uncompress the indexed version of sequences to `$UNIREF_LOCATION`: [uniref90.fasta.dmnd.tar.gz](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90.fasta.dmnd.tar.gz). Or run the following command to download the indexed sequence file into `$UNIREF_LOCATION`:
	`$ metawibele_download_database --database uniref --build uniref90_diamond --install-location $UNIREF_LOCATION`
		* If you are using different version of Diamond, 
			* download raw sequences in fasta format: [uniref90.fasta.tar.gz](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90.fasta.tar.gz). Or run the following command to download the sequence file into `$UNIREF_LOCATION`:
	`$ metawibele_download_database --database uniref --build uniref90_fasta --install-location $UNIREF_LOCATION` 
			* index the sequences using your local Diamond: 
			`$ diamond makedb --in $UNIREF_LOCATION/uniref90.fasta -d $UNIREF_LOCATION/uniref90.fasta`
	* UniRef90 annotation files (5.3 GB): [uniref90_annotations.tar.gz](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90_annotations.tar.gz). Or run the following command to download the annotation files into `$UNIREF_LOCATION`:
	`$ metawibele_download_database --database uniref --build uniref90_annotation --install-location $UNIREF_LOCATION`
		* Annotations are available for the UniRef90 gene families to the following systems:
			* UniProt ID corresponding to the UniRef representative
			* Protein names of the UniRef representative
			* Gene names of the UniRef representative
			* Taxon name and taxon ID of the latest common ancestor (LCA) for each uniref90 cluster
			* Taxon name and taxon ID of the representative of each uniref90 cluster
			* Gene Ontology (GO)
			* KEGG Orthogroups (KOs)
			* EggNOG (including COGs)
			* Pfam domains（Pfams)
		* In most cases, mappings are directly inferred from the annotation of the corresponding UniRef representative sequence in UniProt.

Option 2: Create local uniref databases

* You can also create these databases locally by using MetaWIBELE utility scripts based on the latest release version of UniProt/UniRef, and provide `$DATABASE_LOCATION` as the location to install the database.

* Download and obtain UniProt annotations:

	`$ metawibele_prepare_uniprot_annotation --output $DATABASE_LOCATION`

* Download UniRef sequences and obtain annotations:
	
	`$ metawibele_prepare_uniref_annotation -t uniref90 --output $DATABASE_LOCATION`

* Use `diamond` to index sequences
	
	`$ diamond makedb --in uniref90.fasta -d uniref90.fasta`


##### Domain database
De default, the dependent domain databases **have already been automatically installed** when you install the MetaWIBELE package and you can skip this step. Alternatively, you can also create these domain databases locally and provide `$DATABASE_LOCATION` as the location to install the database.

Create local domain databases (**optional**)

* Download and obtain dependent protein domains information:
	
	`$ metawibele_prepare_domain_databases -t Pfam33.0 --output $DATABASE_LOCATION`
	
	* Pfam domains from the specified version of [Pfam](https://pfam.xfam.org/) database will be downloaded.
	* PDB information from the latest version of [SIFT](https://www.ebi.ac.uk/pdbe/docs/sifts/) database will be downloaded.
	* Domain-domain interactions from version 2.0 of [DOMINE](https://manticore.niehs.nih.gov/cgi-bin/Domine?page=start) database will be downloaded.

	
#### Prepare configuration files
##### Prepare global configuration file
To run MetaWIBELE, you are **required** to customize the global configuration file `metawibele.cfg` and make sure that it's in your working directory. 

* Download `metawibele.cfg` into your working directory:
	* Option 1) run this command to download global configuration file:
		* `$ metawibele_download_config --config-type global`
	* Option 2) obtain copies by right-clicking the link and selecting "save link as":
		* [metawibele.cfg](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/configs/metawibele.cfg)

* Customize your configurations in `metawibele.cfg` before running MetaWIBELE:
	
	Setting the path of dependent databases is the most important to pay attention to and is required for customizing, and most of other sections can leave as defaults.
	
	* Customize the path of dependent databases (**required**):
	
	```
	[database]
	# The absolute path of uniref databases folder.
	uniref_db = 
	# The domain databases used by MetaWIBELE. [data_path] provide the absolute path of the domain databases folder; [none] use the default domain databases installed in the metawibele package [ Default: none ]
	domain_db = none
	```

	* Customize basic information for outputs (e.g. prefix name of output files, etc.) (**optional**):

	```
	[basic]
	# Study name [ Default: MGX ]
	study = MGX
	# The prefix name for output results [ Default: metawibele ]
	basename = metawibele
	``` 
		
	* Customize applied computational resources (e.g. CPU cores, memory, etc.) (**optional**):

	```
	[computation]
	# The number of cores that you’re requesting [ Default: 1 ]
	threads = 1
	# The amount of memory (in MB) that you will be using for your job [ Default: 20000 ] 
	memory = 20000
	# The amount of time (in minute) that you will be using for your job [ Default: 60 ]
	time = 60
	```

	* Customize parameter settings for abundance-based and domain/motif-based annotations (**optional**):

	```
	[abundance]
	# The absolute path of the config file used by MSPminer. [config_file] provide the mspminer config file; [none] use the default config files installed in the metawibele package [ Default: none ]
	mspminer = none
	# The method for normalization [Choices: cpm, relab]. [cpm] copies per million units (sum to 1 million); [relab] relative abundance (sum to 1) [ Default: cpm ]  
	normalize = cpm
	# The minimum abundance for each feature [ Default: 0 ]   
	abundance_detection_level = 0

	[interproscan]
	# Interproscan executable file, e.g. /my/path/interproscan/interproscan.sh [ Default: interproscan.sh ]
	interproscan_cmmd = interproscan.sh
	# The appls used by interproscan: [appls] comma separated list of analyses, [ Choices: CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM,Phobius,SignalP,TMHMM ]; [none] use all all analyses for running [ Default: Pfam,Phobius,SignalP,TMHMM ]
	interproscan_appl = "Pfam,Phobius,SignalP,TMHMM"
	# The number of splitting files which can be annotated in parallel 	[ Default: 1 ]
	split_number = 1
	```
	
	* Customize parameter settings for association with environmental/host phenotypes (**optional**; but you **need** to specify your settings for maaslin2 if you run MetaWIBELE for supervised prioritization, and at least the main phenotype metadata used for prioritization is **required** to set):

	```
	[maaslin2]
	# The absolute path of Maaslin2 executable file, e.g. /my/path/Maaslin2/R/Maaslin2.R [ Default: Maaslin2.R ]
	maaslin2_cmmd = Maaslin2.R
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
	# Apply z-score so continuous metadata are on the same scale [ Default: TRUE ] apply z-score so continuous metadata are on the same scale [ Default: TRUE ]
	standardize = TRUE
	# Generate a heatmap for the significant associations [ Default: FALSE ]
	plot_heatmap = FALSE
	# In heatmap, plot top N features with significant associations [ Default: FALSE ]
	heatmap_first_n = FALSE
	# Generate scatter plots for the significant associations [ Default: FALSE ]
	plot_scatter = FALSE
	# The number of R processes to run in parallel [ Default: 1 ]
	maaslin2_cores = 1
	# The minimum percent of case-control samples used for comparison in which a feature is detected [ Default: 0.1 ]
	tshld_prevalence = 0.10
	# The q-value threshold for significance used as DA annotations [ Default: 0.05 ]
	tshld_qvalue = 0.05
	# The statistic used as effect size [ Choices: coef, mean(log) ]. [coef] represents the coefficient from the model; [mean(log)] represents the difference of mean values between case and control conditions. [  Default: mean(log) ]
	effect_size = mean(log)
	# The main phenotype metadata used for prioritization, e.g. metadata1. [ Default: none ]: skip the association with environmental/phenotypic parameters
	phenotype = none
	# Case and control metadata pairs for phenotype metadata variables; use semicolon to separate variables, e.g. "metadata1:case_status1|control_status1;metadata2:case_status2|control_status2,case_status3|control_status2". [Default: none] where metadata values for the main phenotype will be sorted based on alphabet order and the value with the smallest alphabet order will be treated as control status.
	case_control_status = none
	```


##### Prepare local configuration file
By default, MetaWIBELE will perform by using the local configuration files installed in the package. **Optionally**, you can also make your own local configuration files and provide them with optional arguments to MetaWIBELE. For example, the local characterization configuration file can be provided with `--characterization-config $CHRACTERIZE_CONF` where `$CHRACTERIZE_CONF` is the file including characterization configurations.

* Download local configuration template files (e.g. `characterization.cfg`, `prioritization.cfg `) into your working directory:
	* Option 1) run command line to download local configuration files:
		* `$ metawibele_download_config --config-type local`
	* Option 2) obtain copies by right-clicking the link and selecting "save link as":
		* [characterization.cfg](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/configs/characterization.cfg)
		* [prioritization.cfg](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/configs/prioritization.cfg)
		 
* Modify and provide your own configurations:
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
	# Weight value of prevalence to calculate weighted harmonic mean, named as beta parameter[ Default: 0.50 ] 
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
	# Default: select protein families significantly associated with the main phenotype

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
	
##### Prepare vignette configuration file
**Optionally**, MetaWIBELE can accept user-defined vignette functions of interest for further prioritization. You can make your own vignettes configuration files and provide them with an optional argument to MetaWIBELE. For example, the vignette function configuration file can be provided with `--vignette-config $VIGNETTE_FUNC` where `$VIGNETTE_FUNC` is the file including the functions of interest.

* Download local vignettes template file (`vignettes_function.tsv`) into your working directory:
	* Option 1) Run command line to download vignette configuration files:
		* `$ metawibele_download_config --config-type vignette`
	* Option 2) Obtain copies by right-clicking the link and selecting "save link as":
		* [vignette_function.tsv](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/configs/vignette_function.tsv)

* Make your own configurations:
	* `vignettes_function.tsv` is a tab-separated values file.
	* Two required columns: `type` indicates which type of function it is; `annotation` indicates the specific annotations assigned by MetaWIBELE given an annotation type.
	* Other optional columns: `annotation_type` indicates what type of annotation it is in MetaWIBELE; `description` indicates detailed descriptions of the annotation.
	
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

	`$ metawibele preprocess --help`

	This command yields:

	```
	usage: metawibele [-h] [--version] {characterize,prioritize,preprocess}

	MetaWIBELE workflows: A collection of AnADAMA2 workflows

	positional arguments:
  	{characterize,prioritize,preprocess}	workflow to run

	optional arguments:
  	-h, --help            show this help message and exit
 	 --version             show program's version number and exit
	```
	
* All workflows follow the general command format:

	`$ metawibele $WORKFLOW`

* For specific options of workflow, run:

	`$ metawibele $WORKFLOW --help`
	
	For example: `$ metawibele characterize --help`
	
	This command yields:
	
	```
	usage: characterize.py [-h] [--version] [--threads THREADS]
                       [--characterization-config 	CHARACTERIZATION_CONFIG]
                       [--mspminer-config MSPMINER_CONFIG]
                       [--bypass-clustering] [--bypass-global-homology]
                       [--bypass-domain-motif] [--bypass-interproscan]
                       [--bypass-pfamtogo] [--bypass-domine] [--bypass-sifts]
                       [--bypass-expatlas] [--bypass-psortb]
                       [--bypass-abundance] [--split-number SPLIT_NUMBER]
                       [--bypass-integration] [--study STUDY]
                       [--basename BASENAME] --input-sequence INPUT_SEQUENCE
                       --input-count INPUT_COUNT
                       [--input-metadata INPUT_METADATA] [--output OUTPUT]
                       [-i INPUT] [--config CONFIG] [--local-jobs JOBS]
                       [--grid-jobs GRID_JOBS] [--grid GRID]
                       [--grid-partition GRID_PARTITION]
                       [--grid-benchmark {on,off}]
                       [--grid-options GRID_OPTIONS]
                       [--grid-environment GRID_ENVIRONMENT]
                       [--grid-scratch GRID_SCRATCH] [--dry-run]
                       [--skip-nothing] [--quit-early]
                       [--until-task UNTIL_TASK] [--exclude-task EXCLUDE_TASK]
                       [--target TARGET] [--exclude-target EXCLUDE_TARGET]
                       [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

	A workflow for MetaWIBELE characterization
	```
	

* Run **characterization** workflow

	`$ metawibele characterize --input-sequence <file> --input-count <file> --input-metadata <file> --output <path>`

* Run **prioritization** workflow

	`$ metawibele prioritize --input-annotation <file> --input-attribute <file> --output <path>`

* **Parallelization Options**

	When running any workflow you can add the following command line options to make use of existing computing resources:
	* --local-jobs <1> : Run multiple tasks locally in parallel. Provide the max number of tasks to run at once. The default is one task running at a time.
	* --grid-jobs <0> : Run multiple tasks on a grid in parallel. Provide the max number of grid jobs to run at once. The default is zero tasks are submitted to a grid resulting in all tasks running locally.
	* --grid \<slurm> : Set the grid available on your machine. This will default to the grid found on the machine with options of slurm and sge.
	* --grid-partition \<serial_requeue> : Jobs will be submitted to the partition selected. The default partition selected is based on the default grid.

	For additional workflow options, see the [AnADAMA2](https://github.com/biobakery/anadama2) user manual.


### Standard Workflows
#### MetaWIBELE-characterize
* ##### Input files for for characterization
	* protein sequences for non-redundant gene families (Fasta format file), e.g. [demo_genecatalogs.centroid.faa](https://github.com/biobakery/metawibele/tree/master/examples/input/demo_genecatalogs.centroid.faa)
	* reads counts table for non-redundant gene families (TSV format file), e.g. [demo\_genecatalogs_counts.all.tsv](https://github.com/biobakery/metawibele/tree/master/examples/input/demo_genecatalogs_counts.all.tsv)
	* metadata file (TSV format file), e.g. [demo\_mgx_metadata.tsv](https://github.com/biobakery/metawibele/tree/master/examples/input/demo_mgx_metadata.tsv)
	* the global configuration file in your working directory, e.g. [metawibele.cfg](https://github.com/biobakery/metawibele/tree/master/examples/metawibele.cfg)


* ##### MetaWIBELE-characterize workflow
	`$ metawibele characterize --input-sequence $INPUT_SEQUENCE --input-count $INPUT_COUNT --input-metadata $INPUT_METADATA --output $OUTPUT_DIR`
	
	* Make sure the customized configuration file `metawibele.cfg` is in your working directory.
	* The command replaces `$INPUT_SEQUENCE`, `$INPUT_COUNT`, `$INPUT_METADATA` with three input files, `$OUTPUT_DIR` with the path to the folder to write output files. See the section on **parallelization options** to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings to run all modules. These settings will work for most data sets. However, if you need to customize your workflow settings for the characterization workflow to modify the default settings. You can customize which modules you want to run in your own local configuration file.
		* For example, `--characterization-config $myconfig_file` will modify the default settings when running the characterization modules.
	
* ##### Demo run of MetaWIBELE-characterize

	`$ metawibele characterize --input-sequence demo_genefamilies.centroid.faa --input-count demo_genefamilies_counts.all.tsv --input-metadata demo_mgx_metadata.tsv --output $OUTPUT_DIR`

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
	
	* File name: `$OUTPUT_DIR/$BASENAME_proteinfamilies_annotation.tsv`
	* This file includes the main characterization results (TSV format file).
	* `$OUTPUT_DIR` = the output folder
	* `$BASENAME` = the prefix for output files provided by `metawibele.cfg`
	* This file details the annotation of each protein family in the community. Protein families are groups of evolutionarily-related protein-coding sequences that often perform similar functions.
	* MetaWIBELE annotates protein family by combining global-homology similarity, local-homology similarity, and non-homology based methods.
	* The annotations for each protein family coming from multiple information sources, e.g. biochemical annotation, taxonomic annotation, ecological properties and association with environmental parameters or phenotypes, etc.
			
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
	
	* File name: `$OUTPUT_DIR/$BASENAME_proteinfamilies_annotation.attribute.tsv`
	* This file includes the supplementary results for characterization (TSV format file).
	* Each item is supplemental information about the corresponding results. `AID` is the key to connect `$BASENAME_proteinfamilies_annotation.tsv` with `$BASENAME_proteinfamilies_annotation.attribute.tsv`.
	
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
	
	* File name: `$OUTPUT_DIR/$BASENAME_proteinfamilies_annotation.taxonomy.tsv`
	* This file includes the detailed information about taxonomic annotation (TSV format file).
	* `$BASENAME_proteinfamilies_annotation.taxonomy.all.tsv` shows the taxonomic information at each taxonomic level.
	
	
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
	
	* File name: `$OUTPUT_DIR/$BASENAME_proteinfamilies_nrm.tsv`
	* This file includes the normalized abundance of each protein family across samples (TSV format file).
	* Protein family abundance is reported in copies per million (CPM) units, which is "total sum scaling (TSS)"-style normalization: each sample is constrained to sum to 1 million. First, each protein family is normalized to RPK (reads per kilobase) units for gene length normalization; RPK units reflect relative gene (or transcript) copy number in the community. Then, RPK values are further sum-normalized (CPM) to adjust for differences in sequencing depth across samples. Further information can refer to the normalization approach in [HUMAnN](https://github.com/biobakery/humann). 
	
	**5. Clustering information for protein families**
	
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
	
	* File name: `$OUTPUT_DIR/$BASENAME_proteinfamilies.clstr`
	* This file includes the clustering information for protein families, formatted using an extension-fasta style based on the version of CD-hit clustering file (extension-fasta format file).
	
	**6. Sequences of protein families**  
	
	* File name: `$OUTPUT_DIR/$BASENAME_proteinfamilies.centroid.faa`
	* This file is the protein sequence for representatives of protein families (Fasta format file).

	**7. Intermediate output files**
	
	* Clustering results
		* MetaWIBELE clusters all representative sequences of gene families into protein families. 
		* All intermediate results are in the folder `$OUTPUT_DIR/clustering`.
		
	* Global-homology based search results
		* MetaWIBELE queries each sequence in protein families against the UniRef90 database by performing protein-level search.
		* All intermediate results are in the folder `$OUTPUT_DIR/global_homology_annotation`.
	
	* Domain-motif annotation results
		* MetaWIBELE uses the local-homology approach (domain/motif search) to characterize the secondary structures of protein families.
		* All intermediate results are in the folder `$OUTPUT_DIR/domain_motif_annotation`.
	
	* Abundance-based annotation results
		* MetaWIBELE implements non-homology based strategy compromising (i) taxonomic annotation with phylogenetic binning, (ii) abundance profiling for protein families, and (iii) association with environmental parameters or phenotypes based on differential abundance. 
		* All intermediate results are in the folder `$OUTPUT_DIR/abundance_annotation`.
	

#### MetaWIBELE-prioritize
* ##### Input files for prioritization
	* annotation file produced by MetaWIBELE-characterize workflow (TSV format file):`$BASENAME_proteinfamilies_annotation.tsv`
	* annotation attribute file produced by MetaWIBELE-characterize workflow (TSV format file): `$BASENAME_proteinfamilies_annotation.attribute.tsv`

* ##### MetaWIBELE-prioritize workflow

	`$ metawibele prioritize --input-annotation $INPUT_ANNOTATION --input-attribute $INPUT_ATTRIBUTE --output $OUTPUT_DIR`
	
	* In the command replaces `$INPUT_ANNOTATION`, `$INPUT_ATTRIBUTE` with annotation and attribute files produced by MetaWIBELE-characterize workflow, `$OUTPUT_DIR` with the path to the folder to write output files. 
	* The workflow runs with the default settings for all main tool subtasks. These settings will work for most data sets. However, if you need to customize your workflow settings for the prioritization workflow to determine the optimum setting. Then apply these settings by using options for each task. You can customize your own configuration file.
		* For example, `--prioritization-config $myconfig_file` will modify the default settings when running the prioritization tasks.
	
* ##### Demo run of MetaWIBELE-prioritize

	`$ metawibele prioritize --input-annotation $BASENAME_proteinfamilies_annotation.tsv --input-attribute $BASENAME_proteinfamilies_annotation.attribute.tsv --output $OUTPUT_DIR`

* ##### Output files of MetaWIBELE-prioritize
	**1. unsupervised prioritization**
	
	```
	TID familyID    evidence    value   rank    description note 
	1   Cluster_2   DNA_abundance   6728.564309677419   0.9994472084024323  ranking based on single evidence    
	2   Cluster_2   DNA_prevalence  0.9741935483870968  0.9994472084024323  ranking based on single evidence    
	3   Cluster_2   priority_score  0.9994472084024324  0.9994472084024324  meta ranking based on multiple evidences    
	4   Cluster_269 DNA_abundance   4748.32714451613    0.9983416252072969  ranking based on single evidence    
	5   Cluster_269 DNA_prevalence  0.9483870967741935  0.9964068546158098  ranking based on single evidence    
	6   Cluster_269 priority_score  0.9973733016134971  0.9973733016134971  meta ranking based on multiple evidences 
	...
	```
	
	* File name: 
	`$OUTPUT_DIR/$BASENAME_unsupervised_prioritization.rank.table.tsv`
	* `$OUTPUT_DIR` = the output folder
	* `$BASENAME` = the prefix for output files provided by `metawibele.cfg`
	* This file includes the results of unsupervised prioritization based on ecological properties. Each protein family has a numeric priority score based on meta ranking (TSV format file).
	* `$BASENAME_unsupervised_prioritization.rank.tsv` is the overall ranking for all protein families.

	
	**2. supervised prioritization: numeric ranking**
	
	```
	TID familyID    evidence    value   rank    description note 
	1   Cluster_1058    DNA_within_phenotype_abundance  1455.2607352941177  0.857379767827529   ranking based on single evidence    CD_vs_nonIBD
	2   Cluster_1058    DNA_within_phenotype_prevalence 1.0 1.0 ranking based on single evidence    CD_vs_nonIBD
	3   Cluster_1058    MaAsLin2_DA__mean_log   -3.7441751012074533 0.9336650082918739  ranking based on single evidence    CD_vs_nonIBD
	4   Cluster_1058    MaAsLin2_DA__qvalue 5.01463564253508e-05    0.85437430786268    ranking based on single evidence    CD_vs_nonIBD
	5   Cluster_1058    priority_score  0.9074740723963038  0.9074740723963038  meta ranking based on multiple evidences    CD_vs_nonIBD
	6   Cluster_1152    DNA_within_phenotype_abundance  1617.4002382352937  0.9364289662797125  ranking based on single evidence    CD_vs_nonIBD
	7   Cluster_1152    DNA_within_phenotype_prevalence 1.0 1.0 ranking based on single evidence    CD_vs_nonIBD
	8   Cluster_1152    MaAsLin2_DA__mean_log   -3.7580207457514616 0.9380873410724156  ranking based on single evidence    CD_vs_nonIBD
	9   Cluster_1152    MaAsLin2_DA__qvalue 0.000106062753712042    0.7574750830564784  ranking based on single evidence    CD_vs_nonIBD
	10  Cluster_1152    priority_score  0.8980568683016749  0.8980568683016749  meta ranking based on multiple evidences    CD_vs_nonIBD
	...
	```
	
	* File name: 
		`$OUTPUT_DIR/$BASENAME_supervised_prioritization.rank.table.tsv`
	* * This file includes the results of supervised prioritization by combining ecological properties and environmental/phenotypic properties. Each protein family has a numeric priority score based on meta ranking (TSV format file).
	* `$BASENAME_supervised_prioritization.rank.tsv` is the overall ranking for all protein families.


	**3. supervised prioritization: binary filtering**
	
	***3.1 Select interested subset annotated with at least one of specific biochemical annotations***
	 
	Setting local `$PRIORITIZE_CONF` file as following:
	
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
		* `$ metawibele prioritize --prioritization-config $PRIORITIZE_CONF --bypass-mandatory --selected-output demo_prioritized.selected.tsv --input-annotation $INPUT_ANNOTATION --input-attribute $INPUT_ATTRIBUTE --output $OUTPUT_DIR`
	* Output file name: `$OUTPUT_DIR/demo_prioritized.selected.tsv`
	* This file is the results of binary filtering of protein families based on biochemical annotations (TSV format file).
	* These settings require that each of prioritized protein family should 1) be annotated to domain-domain interaction with the host, and 2) have at least one of the following features: signaling, extracellular, cellWall, outerMembrane, transmembrane 
	
    ***3.2. Select interested subset annotated with multiple specific biochemical annotations simultaneously***
	 
	 Setting `$PRIORITIZE_CONF` as following:
	
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
		* `$ metawibele prioritize --prioritization-config $PRIORITIZE_CONF --bypass-mandatory --selected-output demo_prioritized.selected.tsv --input-annotation $INPUT_ANNOTATION --input-attribute $INPUT_ATTRIBUTE --output $OUTPUT_DIR`
	 * Output file name: `$OUTPUT_DIR/demo_prioritized.selected.tsv`
	 * This file is the results of binary filtering of protein families based on biochemical annotations (TSV format file).
	 * These settings require that each of prioritized protein family should 1) significantly associated with the main phenotype, 2) be annotated to domain-domain interaction with the host, 3) predicted as signal peptides, and 4) have at least one of the following features: extracellular, cellWall, outerMembrane, transmembrane 
	
	***3.3 Select interested subset based on specific functions***
	
	Setting `$PRIORITIZE_CONF` as following:
	
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
		* `$ metawibele prioritize --prioritization-config $PRIORITIZE_CONF --bypass-mandatory --vignette-config my_vignette_function_file  --selected-output demo_prioritized_pilin.tsv --input-annotation $INPUT_ANNOTATION --input-attribute $INPUT_ATTRIBUTE --output $OUTPUT_DIR`
	* Provide your own vignette function file for filtering specific functions.
	* Output file name: `$OUTPUT_DIR/demo_prioritized_pilin.tsv`
	* This file is the result of binary filtering of protein families based on pilin related functions (TSV format file).

***


## Guides to MetaWIBELE Utilities

### Preprocessing sequencing reads to build gene families
A utility workflow in MetaWIBELE package for preprocessing metagenomes reads, used for (i) metagenomic assembly, (ii) gene calling, (iii) gene families (non-redundant gene catalogs) construction, and (iv) gene abundance estimation.

#### Preprocessing workflow
`$ metawibele preprocess --help`

This command yields:

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

A workflow to preprocess shotgun sequencing reads of metagenomes with tasks of metagenomic assembly, gene calling, building gene families, and generating gene abundance for each sample.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --threads THREADS     number of threads/cores for each task to use
  --extension-paired EXTENSION_PAIRED
                        provide the extension for paired fastq files using comma to separate, e.g. .R1.fastq.gz,.R2.fastq.gz | .R1.fastq,.R2.fastq
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
* `--extension-paired` indicates the extension for paired fastq files using comma to separate. It should be specified as ".R1.fastq.gz,.R2.fastq.gz" if the paired fastq files are `$SAMPLE.R1.fastq.gz` and `$SAMPLE.R2.fastq.gz`  
* `--extension` indicates the extension for all fastq files. It should be specified as ".fastq.gz" if the fastq files are `$SAMPLE.fastq.gz` 
* `--output`: the output directory. 

#### Input files of preprocessing workflow
* QC'ed shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format), e.g. 
	- [sample1_R1.fastq.gz](https://github.com/biobakery/metawibele/tree/master/examples/raw_reads/sample1_R1.fastq.gz)
	- [sample1_R2.fastq.gz](https://github.com/biobakery/metawibele/tree/master/examples/raw_reads/sample1_R2.fastq.gz)
	- [sample2_R1.fastq.gz](https://github.com/biobakery/metawibele/tree/master/examples/raw_reads/sample2_R1.fastq.gz)
	- [sample2_R2.fastq.gz](https://github.com/biobakery/metawibele/tree/master/examples/raw_reads/sample2_R2.fastq.gz)
* See the section on parallelization options to optimize the workflow run based on your computing resources. 
* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to modify the default settings, you can change the parameter settings.
	* For example, `--extension-paired "$R1_suffix,$R2_suffix"`, `--extension "$fastq_suffix"` (what are the following part after `$SAMPLE` in the input file names) will modify the default settings when running the assembly task.

#### Demo run of preprocessing workflow

`$ metawibele preprocess --input examples/raw_reads/ --output $OUTPUT_DIR/ --output-basename demo --extension-paired "_R1.fastq.gz,_R2.fastq.gz" --extension ".fastq.gz"`

#### Output files of preprocessing workflow
**Main output files**

The following are the two main output files of the preprocessing utility that are used for MetaWIBELE:

```
$OUTPUT_DIR/$BASENMAE_genecatalogs.centroid.faa
$OUTPUT_DIR/$BASENMAE_genecatalogs_counts.all.tsv

```

**1. demo\_genecatalogs.centroid.faa**
		
```
>PRISM_7938_21108
MKTRRKKQTKRVLAGTLAALMTVSAVPVSNSVVHAEESQDRSELKLRYSSAAPDSYAGWEKWSLPIGNSGIGASVFGGVQ
TERIQLNEKSLWSGGPSDSRPNYNGGNLEEKGKNGQTVKEIQQLFANGDNDAASSKCGELVGLSDDAGVNGYGYYLSYGN
MYLDFKDISDKDVENYERTLDLNTAIAGVEYDNGDTHYTRENFVSYPDNVLVTRLTAEGGDKLNLDVRVEPDNKKGNGSN
NPQPQSYEREWTTNVEDALISIDGQLKDNQMKFSSQTKVLTEGGTTEDGDEKVTVKDAKAVTIITSIGTDYKNDYPVYRT
GESQEQVASRVRAYVDKAADTVEKDSYDTLRQTHVDDYSSIFGRVNLDLGQVPSEKTTDKLLKAYNDGSASDQERRYLEV
...
```
	
* This file provides the amino acid sequences for each gene identified as the representative of its gene family (Fasta format file).
* Each gene is given a MetaWIBELE-specific ID (i.e. PRISM\_7938\_21108) based on the sample and order in which it was identified.


**2. demo\_genecatalogs\_counts.all.tsv**

```
ID	PRISM_7122	PRISM_7147	PRISM_7150	PRISM_7153	PRISM_7184	PRISM_7238	PRISM_7406	PRISM_7408	PRISM_7421
PRISM_7122_03545	42	2	0	4	2	22	1098	117	16
PRISM_7122_03875	197	16	2	15	87	0	0	80	92
PRISM_7122_12067	216	20	6	28	2	12	2006	258	17
PRISM_7122_131770	36	6	6	2	0	70	2274	24	22
PRISM_7122_19039	6	0	0	0	0	14	750	57	5
PRISM_7122_26201	17	6	3	8	2	59	2584	199	20
PRISM_7122_32823	10	0	1	2	0	0	136	19	2
PRISM_7122_38863	16	3	2	8	2	34	905	62	12
PRISM_7122_50124	26	3	3	16	6	35	2139	191	22
...
```

* This file provides the counts of the number of reads mapped to each gene family (rows) in each sample (columns) (TSV format file). 



**Intermediate files** 

**1. assembly results**
	
* `$OUTPUT_DIR/$BASENMAE_contig_sequence.fasta`: contig sequences (Fasta format file).
* The assembly outputs for each sample are in the `$OUTPUT_DIR/assembly/` folder.
	
**2. gene-calling results**
	
* `$OUTPUT_DIR/$BASENMAE_gene_info.tsv`: all gene calls information (TSV format file).
* `$OUTPUT_DIR/$BASENMAE_combined_gene_protein_coding.complete.sorted.fna`: nucleotide sequences for all complete ORFs sorted by gene length (Fasta format file).
* `$OUTPUT_DIR/$BASENMAE_combined_protein.complete.sorted.faa`: protein sequences for all complete ORFs sorted by protein length (Fasta format file).
* `$OUTPUT_DIR/$BASENMAE_combined_gene_protein_coding.sorted.fna`: nucleotide sequences for all ORFs (including partial genes) sorted by gene length (Fasta format file).
* `$OUTPUT_DIR/$BASENMAE_combined_protein.sorted.faa`: protein sequences for all ORFs (including partial genes) sorted by protein length.
* The gene-calling outputs from prodigal are in the `$OUTPUT_DIR/gene_calls` folder. 
* The gene-annotation outputs from prokka are in the `$OUTPUT_DIR/gene_annotation` folder.
	
**3. gene families**
	
* `$OUTPUT_DIR/$BASENMAE_genecatalogs.clstr`: clustering information for non-redundant gene families (extension-fasta format file).
* `$OUTPUT_DIR/$BASENMAE_genecatalogs_counts.all.tsv`: reads counts of gene families per sample (TSV format file).
* `$OUTPUT_DIR/$BASENMAE_genecatalogs.centroid.fna`: nucleotide sequences of representatives for gene families (Fasta format file).
* `$OUTPUT_DIR/$BASENMAE_genecatalogs.centroid.faa`: protein sequences of representatives for gene families.
* All mapping outputs for each sample are in the `$OUTPUT_DIR/mapping` folder. 
	
----

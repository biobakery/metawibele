# MetaWIBELE User Manual

**MetaWIBELE** (**W**orkflow to **I**dentify novel **B**ioactive **Ele**ments in the microbiome) is a workflow to efficiently and systematically identify and prioritize potentially bioactive (and often uncharacterized) gene products in microbial communities. It prioritizes candidate gene products from assembled metagenomes using a combination of sequence homology, secondary-structure-based functional annotations, phylogenetic binning, ecological distribution, and association with environmental parameters or phenotypes to target candidate bioactives.


## Citing MetaWIBELE

**If you use the MetaWIBELE software, please cite our manuscript:**

Yancong Zhang, Amrisha Bhosle, Sena Bae, Lauren J. McIver, Gleb Pishchany, Emma K. Accorsi, Kelsey N. Thompson, Cesar Arze, Ya Wang, Ayshwarya Subramanian, Sean M. Kearney, April Pawluk, Damian R. Plichita, Gholamali Rahnavard, Afrah Shafquat, Ramnik J. Xavier, Hera Vlamakis, Wendy S. Garrett, Andy Krueger, Curtis Huttenhower\*, Eric A. Franzosa\*. "[Discovery of Bioactive Microbial Gene Products in Inflammatory Bowel Disease](https://www.nature.com/articles/s41586-022-04648-7)." *Nature*, 606: 754–760 (2022)

**And feel free to link to MetaWIBELE in your Methods:**

[http://huttenhower.sph.harvard.edu/metawibele](http://huttenhower.sph.harvard.edu/metawibele)


**For additional information, read the [MetaWIBELE Tutorial](https://github.com/biobakery/biobakery/wiki/metawibele)**


If you have questions, please direct it to [the MetaWIBELE channel](https://forum.biobakery.org/c/Microbial-community-profiling/metawibele) of the bioBakery Support Forum.

***

## Contents ##

* [Workflow](#workflow)
	* [Workflow by bypass mode](#workflow-by-bypass-mode) 
* [Install MetaWIBELE](#install-metawibele)
    * [Requirements](#requirements)
    * [Installation](#installation)
    	* [Install MetaWIBELE](#install-metawibele)
    	* [Install databases](#install-databases)
    		* [UniRef database](#uniref-database)
    		* [Domain database](#domain-database)
    	* [Check install](#check-install)
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

* --bypass-clustering
	* do not build protein families
* --bypass-global-homology
	* do not annotate protein families based on global homology information
* --bypass-domain-motif
	* do not annotate protein families based on Interproscan with Phobius/SignalP/TMHMM, Pfam2GO, DOMINE, SIFTS, Expression Atlas database and PSORTb for domain/motif information
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
* --bypass-mspminer
	* do not annotate protein families based on MSPminer information
* --bypass-maaslin
	* do not annotate protein families based on MaAsLin2 information
* --bypass-integration
	* do not integrate annotations for protein families
* --bypass-mandatory
	* do not prioritize protein families to calculate continuous priority scores
* --bypass-optional
	* do not prioritize protein families based on selecting our for interested annotations (optional prioritization)


## Install MetaWIBELE
### Requirements
1. [Python](https://www.python.org/) (version >= 3.6, requiring numpy, pandas, scipy packages; *tested 3.6, 3.7*)
2. [AnADAMA2](https://huttenhower.sph.harvard.edu/anadama2) (version >= 0.7.4; *tested 0.7.4, 0.8.0*)
3. [CD-hit](http://weizhongli-lab.org/cd-hit/) (version >= 4.7; *tested 4.7*)
4. [Diamond](http://www.diamondsearch.org/index.php) (version >= 0.9.24; *tested 0.9.24*)
5. [MSPminer](https://academic.oup.com/bioinformatics/article/35/9/1544/5106712?login=true) (version >= 1.0.0; licensed software; *tested version 1.0.0*)
6. [MaAsLin2](https://huttenhower.sph.harvard.edu/maaslin2) (version >= 1.5.1; *tested 1.5.1*)
7. [Interproscan](https://interproscan-docs.readthedocs.io/en/latest/) (version >= 5.31-70) (installing with activating Phobius/SignalP/TMHMM analyses; InterProScan 5.51-85.0 or later are recommended for potential simpler installation; *tested 5.31-70, 5.51-85.0*; you can skip to install Interproscan if you’d like to ignore domain/motif annotation by running MetaWIBELE with “--bypass-interproscan”)
8. [Signalp](http://www.cbs.dtu.dk/services/SignalP-4.1/) (version >= 4.1; licensed software; installing with integrating in interproscan; used for domain/motif annotation; *tested 4.1*)
9. [TMHMM](http://www.cbs.dtu.dk/services/TMHMM/) (version >= 2.0; licensed software; installing with integrating in interproscan; used for domain/motif annotation; *tested 2.0*)
10. [Phobius](http://phobius.sbc.su.se/) (version >= 1.01; licensed software; installing with integrating in interproscan; used for domain/motif annotation; *tested 1.01*)
11. [PSORTb](https://psort.org/documentation/index.html) (version >= 3.0) (licensed software; used for domain/motif annotation; *tested 3.0*; you can skip to install psortb if you’d like to ignore protein subcellular localization annotations by running MetaWIBELE with “--bypass-psortb”)
12. **Optional**: only required if using MetaWIBELE utilities to prepare inputs for MetaWIBELE using metagenomic sequencing reads
	* [MEGAHIT](https://github.com/voutcn/megahit) (version >= 1.1.3; tested *1.1.3*) 
	* [Prokka](https://github.com/tseemann/prokka) (version >= 1.14-dev; recommend to not set '-c' parameter when running prodigal with metagenome mode; *tested 1.14-dev*)
	* [Prodigal](https://github.com/hyattpd/Prodigal) (version >= 2.6; *tested 2.6*)
	* [USEARCH](http://www.drive5.com/usearch/) (version >= 9.0.2132; licensed software; *tested 9.0.2132*)
	* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.3.2; *tested 2.3.2*)
	* [SAMtools](https://github.com/samtools/) (version >= 1.9; *tested 1.9*)
	* [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) (version >= 1.6.2; *tested 1.6.2*)

**Note:** Please install the required software in a location in your `$PATH`. If you always run with gene families (non-redundant gene catalogs), the optional softwares are not required. Also if you always run with one or more bypass options (for information on bypass options, see optional arguments to the section [Workflow by bypass mode](#workflow-by-bypass-mode)). The software required for the steps you bypass does not need to be installed.


### Installation
#### Install MetaWIBELE
You only need to do **any one** of the following options to install the MetaWIBELE package.

**Option 1: Installing with docker**

* `$ docker pull biobakery/metawibele`
* This docker image includes most of the dependent software packages.
* Large software packages and those with licenses are **NOT** included in this image and needed to be installed additionally:
	* Users should review the license terms and install these packages manually. 
	* Softwares with the license : [MSPminer](https://academic.oup.com/bioinformatics/article/35/9/1544/5106712?login=true), [Signalp](http://www.cbs.dtu.dk/services/SignalP-4.1/), [TMHMM](http://www.cbs.dtu.dk/services/TMHMM/), [Phobius](http://phobius.sbc.su.se/), [PSORTb](https://psort.org/documentation/index.html)
	* Softwares with large size: [Interproscan](https://interproscan-docs.readthedocs.io/en/latest/) (**Note:** We recommen dinstalling InterProScan 5.51-85.0 (requiring at least Java 11) or later for potential simpler installation, and active Phobius/SignalP/TMHMM analyses by customizing your `interproscan.properties` configuration, see more details from [InterProScan document](https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html)).
	

**Option 2: Installing with pip**

* `$ pip install metawibele`
* If you do not have write permissions to `/usr/lib/`, then add the option --user to the install command. This will install the python package into subdirectories of `~/.local/`. Please note when using the --user install option on some platforms, you might need to add `~/.local/bin/` to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `metawibele: command not found` when trying to run MetaWIBELE after installing with the --user option.
	
**Option 3: Installing with conda**

* `$ conda install -c biobakery metawibele`
 

#### Install databases
To run metawibele, you need to install the dependent databases: 1) uniref databases (**required**); 2) domain databases (optional).

##### 1. UniRef database
UniRef database is **required** if you will use MetaWIBELE to do global-homology based annotation and taxonomic annotation. You can use **any one** of the following options to install the database with providing `$UNIREF_LOCATION` as the location for installation.

**NOTE:** Please point to this location for the default uniref database in the global config file (`metawibele.cfg`, see details in the section [Prepare global configuration file](#prepare-global-configuration-file)). Alternatively, you can either set the location with the environment variable `$UNIREF_LOCATION`, or move the downloaded database to the folder named `uniref_database` in the current working directory.

**Option 1: Download uniref databases (Recommended)**

We have built the dependent UniRef database based on UniProt/UniRef 2019_01 sequences and annotations. You can download and uncompress this database (both sequences and annotations) and provide `$UNIREF_LOCATION` as the location to install the database.
	
* UniRef90 sequence file (20 GB): 	
	* If you are using Diamond v0.9.24, just download and uncompress the indexed version of sequences to `$UNIREF_LOCATION`: [uniref90.fasta.dmnd.tar.gz](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90.fasta.dmnd.tar.gz) 
	* Alternatively, run the following command to download the indexed sequence file by Diamond v0.9.24 into `$UNIREF_LOCATION`:
	`$ metawibele_download_database --database uniref --build uniref90_diamond --install-location $UNIREF_LOCATION`
	* If you are using different version of Diamond, download raw sequences in fasta format: [uniref90.fasta.tar.gz](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90.fasta.tar.gz) 
		* Alternatively, run the following command to download the sequence file into `$UNIREF_LOCATION`:
	`$ metawibele_download_database --database uniref --build uniref90_fasta --install-location $UNIREF_LOCATION` 
		* And then, index the sequences using your local Diamond: 
			`$ diamond makedb --in $UNIREF_LOCATION/uniref90.fasta -d $UNIREF_LOCATION/uniref90.fasta`

* UniRef90 annotation files (2.8 GB): 
	* Download the annotation files and uncompress into `$UNIREF_LOCATION`: [uniref90_annotations.tar.gz](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/uniref90_annotations.tar.gz). 
	* Alternatively, run the following command to download the annotation files into `$UNIREF_LOCATION`:
	`$ metawibele_download_database --database uniref --build uniref90_annotation --install-location $UNIREF_LOCATION`


**Option 2: Create local uniref databases**

* You can also create these databases locally by using MetaWIBELE utility scripts based on the latest release version of UniProt/UniRef, and provide `$UNIREF_LOCATION` as the location to install the database.

* Download and obtain UniProt annotations:

	`$ metawibele_prepare_uniprot_annotation --output $UNIREF_LOCATION `

* Download UniRef sequences and obtain annotations:
	
	`$ metawibele_prepare_uniref_annotation -t uniref90 --output $UNIREF_LOCATION `

* Use `diamond` to index sequences
	
	`$ diamond makedb --in uniref90.fasta -d uniref90.fasta`



##### 2. Domain database
De default, the dependent domain databases **have already been automatically installed** when you install the MetaWIBELE package and you can skip this step. Alternatively, you can also create these domain databases locally and provide `$DOMAIN_LOCATION` as the location to install the database.

**NOTE:** Please point to this location for the default domain database in the global config file (`metawibele.cfg`, see details in the section [Prepare global configuration file](#prepare-global-configuration-file)).

Create local domain databases (**optional**)

* Download and obtain dependent protein domains information:
	
	`$ metawibele_prepare_domain_databases -t Pfam33.0 --output $DOMAIN_LOCATION `
	
	* Pfam domains from the specified version of [Pfam](https://pfam.xfam.org/) database will be downloaded.
	* PDB information from the latest version of [SIFT](https://www.ebi.ac.uk/pdbe/docs/sifts/) database will be downloaded.
	* Domain-domain interactions from version 2.0 of [DOMINE](https://manticore.niehs.nih.gov/cgi-bin/Domine?page=start) database will be downloaded.


#### Check install
To check out the install of MetaWIBELE packages and all dependencies (tools and databases), run the command: 
	
* `$ metawibele_check_install`
* By default, it will check both MetaWIBEKE packages and all dependencies.
* Alternatively, add the option "--types {metawibele,required,optional,all}" to check the install of specific package or dependency.
 
	
#### Prepare configuration files
##### Prepare global configuration file
To run MetaWIBELE, you are **required** to customize the global configuration file `metawibele.cfg` and make sure that it's in the current working directory. 

**NOTE:** De default, MetaWIBELE will use the global configurations from `metawibele.cfg` in the current working directory. Alternatively you can always provide the location of the global configuration file you would like to use with the "--global-config " option to metawibele (see more in the section [How to run](#how-to-run)).

* Download `metawibele.cfg` into your current working directory by any one of the following options:
	* Option 1) obtain copies by right-clicking the link and selecting "save link as": [metawibele.cfg](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/configs/metawibele.cfg)
	* Option 2) run this command to download global configuration file:
	`$ metawibele_download_config --config-type global`

* Customize your configurations in `metawibele.cfg` before running MetaWIBELE:
	
	The path of the dependent databases is **required** for customization. Most of other sections can be left as defaults.
	
	* Customize the path of dependent databases (**required**):
	
	```
	[database]
	# The absolute path of uniref databases folder.
	uniref_db = 
	# The domain databases used by MetaWIBELE. [data_path] provide the absolute path of the domain databases folder; [none] use the default domain databases installed in the metawibele package. [ Default: none ]
	domain_db = none
	```

	* Customize basic information for outputs (e.g. prefix name of output files, etc.) (**optional**):

	```
	[basic]
	# Study name. [ Default: MGX ]
	study = MGX
	# The prefix name for output results. [ Default: metawibele ]
	basename = metawibele
	``` 
		
	* Customize applied computational resources (e.g. CPU cores, memory, etc.) (**optional**):

	```
	[computation]
	# The number of cores that you’re requesting. [ Default: 1 ]
	threads = 1
	# The amount of memory (in MB) that you will be using for your job. [ Default: 20000 ] 
	memory = 20000
	# The amount of time (in minute) that you will be using for your job. [ Default: 60 ]
	time = 60
	```

	* Customize parameter settings for abundance-based and domain/motif-based annotations (**optional**):

	```
	[abundance]
	# The absolute path of the config file used by MSPminer. [config_file] provide the mspminer config file; [none] use the default config files installed in the metawibele package. [ Default: none ]
	mspminer = none
	# The method for normalization [Choices: cpm, relab]. [cpm] copies per million units (sum to 1 million); [relab] relative abundance (sum to 1). [ Default: cpm ]  
	normalize = cpm
	# The minimum abundance for each feature [ Default: 0 ]   
	abundance_detection_level = 0

	[msp]
	# The minimum fraction of taxonomy classified genes in each MSP [0-1]. [Default: 0.10]
	tshld_classified = 0.10 
	# The minimum percent differences between the most and second dominant taxon for each MSP [0-1]. [Default: 0.50] 
	tshld_diff = 0.50 

	[interproscan]
	# Interproscan executable file, e.g. /my/path/interproscan/interproscan.sh [ Default: interproscan.sh ]
	interproscan_cmmd = interproscan.sh
	# The appls used by interproscan: [appls] comma separated list of analyses, [ Choices: CDD,COILS,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM,Phobius,SignalP,TMHMM ]; [all] use all all analyses for running. [ Default: all ]
	interproscan_appl = all
	# The number of splitting files which can be annotated in parallel 	[ Default: 1 ]
	split_number = 1
	```
	
	* Customize parameter settings for association with environmental/host phenotypes (**optional**; but you are **required** to specify your settings for [MaAsLin2](https://huttenhower.sph.harvard.edu/maaslin2) if you run MetaWIBELE for supervised prioritization, and at least the main `phenotype` metadata used for prioritization is **required** to set):

	```
	[maaslin2]
	# The absolute path of Maaslin2 executable file, e.g. /my/path/Maaslin2/R/Maaslin2.R [ Default: Maaslin2.R ]
	maaslin2_cmmd = Maaslin2.R
	# The minimum abundance for each feature. [ Default: 0 ]  
	min_abundance = 0
	# The minimum percent of samples for which a feature is detected at minimum abundance. [ Default: 0.1 ]
	min_prevalence = 0.1
	# Keep features with variance greater than. [Default: 0.0]
	min_variance = 0
	# The q-value threshold for significance. [ Default: 0.25 ]
	max_significance = 0.25
	# The normalization method to apply. [ Choices: TSS, CLR, CSS, NONE, TMM ], [ Default: TSS ]
	normalization = NONE
	# The transform to apply [ Choices: LOG, LOGIT, AST, NONE ].  [ Default: LOG ]
	transform = LOG
	# The analysis method to apply [ Choices: LM, CPLM, ZICP, NEGBIN, ZINB ]. [ Default: LM ]
	analysis_method = LM
	# The fixed effects for the model, comma-delimited for multiple effects. [ Default: all ]
	fixed_effects = all
	# The random effects for the model, comma-delimited for multiple effects. [ Default: none ]
	random_effects = none
	# The correction method for computing the q-value. [ Default: BH ]
	correction = BH
	# Apply z-score so continuous metadata are on the same scale [ Default: TRUE ] apply z-score so continuous metadata are on the same scale. [ Default: TRUE ]
	standardize = TRUE
	# Generate a heatmap for the significant associations. [ Default: FALSE ]
	plot_heatmap = FALSE
	# In heatmap, plot top N features with significant associations. [ Default: FALSE ]
	heatmap_first_n = FALSE
	# Generate scatter plots for the significant associations. [ Default: FALSE ]
	plot_scatter = FALSE
	# The number of R processes to run in parallel. [ Default: 1 ]
	maaslin2_cores = 1
	# The factor to use as a reference for a variable with more than two levels provided as a string of 'variable,reference' semi-colon delimited for multiple variables. NOTE: A space between the variable and reference will not error but will cause an inaccurate result. [ Default: NA ]
	reference = NA
	# The minimum percent of case-control samples used for comparison in which a feature is detected. [ Default: 0.1 ]
	tshld_prevalence = 0.10
	# The q-value threshold for significance used as DA annotations. [ Default: 0.05 ]
	tshld_qvalue = 0.05
	# The statistic used as effect size [ Choices: coef, mean(log) ]. [coef] represents the coefficient from the model; [mean(log)] represents differences of mean log-scaled abundances between case and control conditions. [  Default: mean(log) ]
	effect_size = mean(log)
	# The main phenotype metadata used for prioritization, e.g. metadata1. [ Default: none ]: skip the association with environmental/phenotypic parameters
	phenotype = none
	```


##### Prepare local configuration file
By default, MetaWIBELE will perform by using the local configuration files installed in the package. **Optionally**, you can also make your own local configuration files and provide them with optional arguments to MetaWIBELE. For example, the local characterization configuration file can be provided with `--characterization-config $CHRACTERIZE_CONF` where `$CHRACTERIZE_CONF` is the file including characterization configurations.

* Download local configuration template files (e.g. `characterization.cfg`, `prioritization.cfg `) into your working directory:
	* Option 1) run command line to download local configuration files:
		* `$ metawibele_download_config --config-type local`
	* Option 2) obtain copies by right-clicking the link and selecting "save link as":
		* [characterization.cfg](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/configs/characterization.cfg)
		* [prioritization.cfg](http://huttenhower.sph.harvard.edu/MetaWIBELE_data/configs/prioritization.cfg)
		 
* Customize local configurations:
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

	`$ metawibele --help`

	This command yields:

	```
	usage: metawibele [-h] [--global-config GLOBAL_CONFIG]
                  {characterize,prioritize,preprocess}

	MetaWIBELE workflows: A collection of AnADAMA2 workflows

	positional arguments:
  	{profile,characterize,prioritize,preprocess}
                        workflow to run

	optional arguments:
  	-h, --help            show this help message and exit
  	--global-config GLOBAL_CONFIG the global configuration file of MetaWIBELE (default: None)
	```
	
**NOTE:** De default, MetaWIBELE will use the global configurations from `metawibele.cfg` in the current working directory. Alternatively, you can always provide the location of the global configuration file you would like to use with the "--global-config " option to metawibele.	
	
* All workflows follow the general command format:

	`$ metawibele $WORKFLOW`

* For specific options of workflow, run:

	`$ metawibele $WORKFLOW --help`
	
	For example: `$ metawibele characterize --help`
	
	This command yields:
	
	```
	usage: characterize.py [-h] [--version] [--threads THREADS]
                       [--characterization-config CHARACTERIZATION_CONFIG]
                       [--mspminer-config MSPMINER_CONFIG]
                       [--bypass-clustering] [--bypass-global-homology]
                       [--bypass-domain-motif] [--bypass-interproscan]
                       [--bypass-pfamtogo] [--bypass-domine] [--bypass-sifts]
                       [--bypass-expatlas] [--bypass-psortb]
                       [--bypass-abundance] [--bypass-mspminer] [--bypass-maaslin]
                       [--split-number SPLIT_NUMBER]
                       [--bypass-integration] [--study STUDY]
                       [--basename BASENAME] --input-sequence INPUT_SEQUENCE
                       --input-count INPUT_COUNT --input-metadata
                       INPUT_METADATA [--output OUTPUT] [-i INPUT]
                       [--config CONFIG] [--local-jobs JOBS]
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
	* protein sequences for non-redundant gene families (Fasta format file), e.g. [demo_genecatalogs.centroid.faa](https://raw.githubusercontent.com/biobakery/metawibele/master/examples/input/demo_genecatalogs.centroid.faa)
	* reads counts table for non-redundant gene families (TSV format file), e.g. [demo\_genecatalogs_counts.all.tsv](https://raw.githubusercontent.com/biobakery/metawibele/master/examples/input/demo_genecatalogs_counts.all.tsv)
	* metadata file (TSV format file), e.g. [demo\_mgx_metadata.tsv](https://raw.githubusercontent.com/biobakery/metawibele/master/examples/input/demo_mgx_metadata.tsv)
	* the global configuration file in the current working directory, e.g. [metawibele.cfg](https://raw.githubusercontent.com/biobakery/metawibele/master/examples/metawibele.cfg)
	

* ##### MetaWIBELE-characterize workflow
	`$ metawibele characterize --input-sequence $INPUT_SEQUENCE --input-count $INPUT_COUNT --input-metadata $INPUT_METADATA --output $OUTPUT_DIR`
	
	* Make sure the customized configuration file `metawibele.cfg` is in your working directory.
	* The command replaces `$INPUT_SEQUENCE`, `$INPUT_COUNT`, `$INPUT_METADATA` with three input files, `$OUTPUT_DIR` with the path to the folder to write output files. See the section on **parallelization options** to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings to run all modules. These settings will work for most data sets. However, if you need to customize your workflow settings for the characterization workflow to modify the default settings. You can customize which modules you want to run in your own local configuration file.
		* For example, `--characterization-config $myconfig_file` will modify the default settings when running the characterization modules.
	
* ##### Demo run of MetaWIBELE-characterize

	`$ metawibele characterize --input-sequence demo_genecatalogs.centroid.faa --input-count demo_genecatalogs_counts.all.tsv --input-metadata demo_mgx_metadata.tsv --output $OUTPUT_DIR`

* ##### Output files of MetaWIBELE-characterize
	**1. Annotation file**
	
	```
	familyID	annotation	feature	category	method	AID
	Cluster_1	good	quality	note	Quality_control	NA
	Cluster_1	demo	study	project	Shotgun	NA
	Cluster_1	Cluster_1	protein_family	Denovo_clustering	CD-hit	Cluster_1__Denovo_clustering
	Cluster_1	UniRef90_A7B522	strong_homology	UniRef90_homology	UniRef90	Cluster_1__UniRef90_homology
	Cluster_1	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	Terminal	Taxonomy_characterization	Taxonomy_annotation	Cluster_1__Taxonomy_characterization
	Cluster_1	4734.313443	DNA-CD_abundance	Denovo_characterization	DNA	Cluster_1__DNA-CD_abundance
	Cluster_1	0.926470588	DNA-CD_prevalence	Denovo_characterization	DNA	Cluster_1__DNA-CD_prevalence
	Cluster_1	1048.199147	DNA-nonIBD_abundance	Denovo_characterization	DNA	Cluster_1__DNA-nonIBD_abundance
	Cluster_1	1	DNA-nonIBD_prevalence	Denovo_characterization	DNA	Cluster_1__DNA-nonIBD_prevalence
	Cluster_1	4734.313443	DNA-within-phenotype_abundance	Denovo_characterization	DNA	Cluster_1__DNA-within-phenotype_abundance
	Cluster_1	1	DNA-within-phenotype_prevalence	Denovo_characterization	DNA	Cluster_1__DNA-within-phenotype_prevalence
	Cluster_1	3505.608677	DNA_abundance	Denovo_characterization	DNA	Cluster_1__DNA_abundance
	Cluster_1	0.950980392	DNA_prevalence	Denovo_characterization	DNA	Cluster_1__DNA_prevalence
	...
	```
	
	* File name: `$OUTPUT_DIR/finalized/$BASENAME_proteinfamilies_annotation.tsv`
	* This file includes the main characterization results (TSV format file).
	* `$OUTPUT_DIR` = the output folder
	* `$BASENAME` = the prefix for output files provided by `metawibele.cfg`
	* This file details the annotation of each protein family in the community. Protein families are groups of evolutionarily-related protein-coding sequences that often perform similar functions.
	* MetaWIBELE annotates protein family by combining global-homology similarity, local-homology similarity, and non-homology based methods.
	* The annotations for each protein family coming from multiple information sources, e.g. biochemical annotation, taxonomic annotation, ecological properties and association with environmental parameters or phenotypes, etc.
			
	**2. Attribute file**
	
	```
	TID	AID	key	value
	1	Cluster_1__Denovo_clustering	repID	PRISM_7861_211956
	2	Cluster_1__Denovo_clustering	rep_length	85
	3	Cluster_1__Denovo_clustering	cluster_size	24
	4	Cluster_1__UniRef90_homology	UniProtKB	A7B522
	5	Cluster_1__UniRef90_homology	Protein_names	Translation initiation factor IF-1
	6	Cluster_1__UniRef90_homology	query_cov_type	high_confidence
	7	Cluster_1__UniRef90_homology	mutual_cov_type	high_confidence
	8	Cluster_1__UniRef90_homology	identity	100
	9	Cluster_1__UniRef90_homology	query_coverage	1
	10	Cluster_1__UniRef90_homology	mutual_coverage	1
	11	Cluster_1__UniRef90_homology	taxa_id	1
	12	Cluster_1__UniRef90_homology	taxa_name	Unclassified
	13	Cluster_1__Taxonomy_characterization	taxa_id	411470
	...
	```
	
	* File name: `$OUTPUT_DIR/finalized/$BASENAME_proteinfamilies_annotation.attribute.tsv`
	* This file includes the supplementary results for characterization (TSV format file).
	* Each item is supplemental information about the corresponding results. `AID` is the key to connect `$BASENAME_proteinfamilies_annotation.tsv` with `$BASENAME_proteinfamilies_annotation.attribute.tsv`.
	
	**3. Taxonomic file**
	
	```
	familyID	study	map_type	query_type	mutual_type	identity	query_coverage	mutual_coverage	detail	Tax	TaxID	Rep_Tax	Rep_TaxID	UniProtKB	unirefID	note	msp_name	msp_taxa_name	msp_taxa_id	MSP_Tax	MSP_TaxID	MSP_Rep_Tax	MSP_Rep_TaxID	taxa_id	taxa_name	taxa_rank	taxa_lineage
	Cluster_1	demo	UniRef90_characterized	high_confidence	high_confidence	100	1	1	Translation initiation factor IF-1	root	1	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	A7B522	UniRef90_A7B522	good	msp_1	Ruminococcus gnavus	33038	Lachnospiraceae	186803	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	411470	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	Terminal	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_gnavus|t__Ruminococcus_gnavus_(strain_ATCC_29149_/_VPI_C7-9)
	Cluster_10	demo	UniRef90_uncharacterized	high_confidence	high_confidence	97.2	0.829457364	0.829457364	Uncharacterized protein	Clostridiales	186802	[Eubacterium] rectale	39491	A0A0M6W8Q5	UniRef90_A0A0M6W8Q5	good	msp_2	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	39491	[Eubacterium] rectale	Species	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|unclassified_Lachnospiraceae|s__[Eubacterium]_rectale
	Cluster_100	demo	UniRef90_characterized	high_confidence	high_confidence	99.7	1	1	S-adenosylmethionine:tRNA ribosyltransferase-isomerase	Clostridiales	186802	[Eubacterium] rectale	39491	A0A0M6WLI9	UniRef90_A0A0M6WLI9	good	msp_2	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	39491	[Eubacterium] rectale	Species	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|unclassified_Lachnospiraceae|s__[Eubacterium]_rectale
	Cluster_1000	demo	UniRef90_uncharacterized	high_confidence	high_confidence	99.5	1	1	Domain of uncharacterized function (DUF1836)	Clostridiales	186802	[Eubacterium] rectale	39491	A0A174GTW6	UniRef90_A0A174GTW6	good	msp_2	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	39491	[Eubacterium] rectale	Species	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|unclassified_Lachnospiraceae|s__[Eubacterium]_rectale
	Cluster_1001	demo	UniRef90_characterized	high_confidence	high_confidence	97.7	1	1	DUF624 domain-containing protein	Clostridiales	186802	[Eubacterium] rectale	39491	A0A173U990	UniRef90_A0A173U990	good	msp_2	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	[Eubacterium] rectale	39491	39491	[Eubacterium] rectale	Species	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|unclassified_Lachnospiraceae|s__[Eubacterium]_rectale
	Cluster_1002	demo	UniRef90_characterized	high_confidence	high_confidence	93.5	1	0.981900452	Proton-coupled thiamine transporter YuaJ	Clostridiales	186802	Ruminococcus gnavus	33038	A0A2N5PZR6	UniRef90_A0A2N5PZR6	good	msp_unknown	NA	NA	Ruminococcus gnavus	33038	Ruminococcus gnavus	33038	33038	Ruminococcus gnavus	Species	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_gnavus
	Cluster_1003	demo	UniRef90_characterized	high_confidence	high_confidence	98.6	1	1	HAD hydrolase, family IA, variant 1	Clostridiales	186802	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	A7B2S9	UniRef90_A7B2S9	good	msp_1	Ruminococcus gnavus	33038	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	411470	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	Terminal	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_gnavus|t__Ruminococcus_gnavus_(strain_ATCC_29149_/_VPI_C7-9)
	Cluster_1004	demo	UniRef90_uncharacterized	high_confidence	high_confidence	99.5	1	0.96	YheO-like protein	Clostridiales	186802	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	A7B3N1	UniRef90_A7B3N1	good	msp_1	Ruminococcus gnavus	33038	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	411470	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	Terminal	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_gnavus|t__Ruminococcus_gnavus_(strain_ATCC_29149_/_VPI_C7-9)
	Cluster_1005	demo	UniRef90_characterized	high_confidence	high_confidence	100	1	1	Zf-HC2 domain-containing protein	Clostridiales	186802	Ruminococcus gnavus	33038	A0A2N5Q1F5	UniRef90_A0A2N5Q1F5	good	msp_1	Ruminococcus gnavus	33038	Ruminococcus gnavus	33038	Ruminococcus gnavus	33038	33038	Ruminococcus gnavus	Species	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_gnavus
	Cluster_1006	demo	UniRef90_uncharacterized	high_confidence	high_confidence	99.5	1	1	YigZ family protein	Clostridiales	186802	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	A7B7M9	UniRef90_A7B7M9	good	msp_1	Ruminococcus gnavus	33038	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	411470	411470	Ruminococcus gnavus (strain ATCC 29149 / VPI C7-9)	Terminal	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_gnavus|t__Ruminococcus_gnavus_(strain_ATCC_29149_/_VPI_C7-9)
	Cluster_1007	demo	UniRef90_characterized	high_confidence	high_confidence	99.5	1	1	Sugar phosphate isomerase/epimerase	Clostridiales	186802	Ruminococcus gnavus	33038	A0A2N5NPC6	UniRef90_A0A2N5NPC6	good	msp_1	Ruminococcus gnavus	33038	Ruminococcus gnavus	33038	Ruminococcus gnavus	33038	33038	Ruminococcus gnavus	Species	k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Ruminococcus_gnavus
	...
	```
	
	* File name: `$OUTPUT_DIR/finalized/$BASENAME_proteinfamilies_annotation.taxonomy.tsv`
	* This file includes the detailed information about taxonomic annotation (TSV format file).
	* `$BASENAME_proteinfamilies_annotation.taxonomy.all.tsv` shows the taxonomic information at each taxonomic level.
	
	
	**4. Abundance file**
	
	```
	ID	PRISM_7122	PRISM_7147	PRISM_7150	PRISM_7153	PRISM_7184	PRISM_7238	PRISM_7406	PRISM_7408	PRISM_7421	PRISM_7445	PRISM_7486	PRISM_7547	PRISM_7658	PRISM_7662	PRISM_7744	PRISM_7759	PRISM_7791	PRISM_7843	PRISM_7847	PRISM_7855	PRISM_7858	PRISM_7860	PRISM_7861	PRISM_7862	PRISM_7870	PRISM_7874	PRISM_7875	PRISM_7879	PRISM_7899	PRISM_7904	PRISM_7906	PRISM_7908	PRISM_7909	PRISM_7910	PRISM_7911	PRISM_7912	PRISM_7938	PRISM_7941	PRISM_7947	PRISM_7948	PRISM_7955	PRISM_7971	PRISM_7989	PRISM_8095	PRISM_8226	PRISM_8244	PRISM_8264	PRISM_8283	PRISM_8332	PRISM_8336	PRISM_8361	PRISM_8374	PRISM_8377	PRISM_8406	PRISM_8452	PRISM_8462	PRISM_8466	PRISM_8467	PRISM_8475	PRISM_8483	PRISM_8485	PRISM_8496	PRISM_8523	PRISM_8534	PRISM_8537	PRISM_8550	PRISM_8564	PRISM_8565	PRISM_8573	PRISM_8577	PRISM_8589	PRISM_8591	PRISM_8592	PRISM_8624	PRISM_8629	PRISM_8675	PRISM_8683	PRISM_8746	PRISM_8749	PRISM_8753	PRISM_8754	PRISM_8758	PRISM_8764	PRISM_8765	PRISM_8774	PRISM_8776	PRISM_8783	PRISM_8784	PRISM_8788	PRISM_8789	PRISM_8794	PRISM_8800	PRISM_8802	PRISM_8806	PRISM_8807	PRISM_8841	PRISM_8843	PRISM_8847	PRISM_8878	PRISM_8892	PRISM_9126	PRISM_9148
	Cluster_1	1628.57	344.191	1779.04	361.815	21274.3	2531.96	1758.21	1325.19	107.936	2479.62	75.0717	10288	140.489	526.223	878.457	544.58	23741.5	23477.1	1020.34	612.181	1677.33	798.845	249.601	350.223	1297.56	565.739	349.645	2893.53	310.445	242.589	1151.99	397.229	1296.43	157.584	1509.21	527.071	730.622	2006.88	1762.46	88.4127	100.735	0	1512.46	1975.95	1348.34	1816.61	926.93	282.481	124.636	0	315.456	1730.5	1576.04	145.977	0	22840.6	276.333	594.827	1821.89	2485.37	11385.1	7361.11	2019.45	1551.36	984.562	16840.9	3328	111359	404.952	424.179	109.963	2260.72	396.53	1677.01	4283.67	1771.56	156.755	0	1171.73	420.518	1235.05	1441.67	4171.22	3716.82	488.319	724.294	1717.69	539.79	529.208	1178.29	144.675	85.8857	2293.57	0	1391.19	1253.93	5509.35	2372.51	1608.11	1249.66	2349.23	2531.25
	Cluster_10	0	1775.32	1406.68	1703.73	1168.16	61.7906	0	201.504	1249.7	0	1723.06	376.605	1961.33	2162	1653.13	837.274	810.55	0	1859.61	1337.51	1698.26	1445.77	1717.21	1444.15	1019.4	2490.13	230.386	1463.19	2277.4	2028.38	2049.48	1046.96	1952.54	1682.12	520.897	1835.7	1752.66	3.03991	331.803	865.525	1766.34	9685.46	0	72.3325	0	1455.8	0	1669.18	1748.08	0	2188.82	23.4298	0	1682.03	0	106.738	1480.39	1383.32	1200.47	0	0	0	0	17.3256	0	0	0	643.647	1067.32	1816.73	1449.12	0	1704.54	0	0	14.4112	2326.57	4110.52	22.3788	1080.07	1348.1	0	1624.1	1662.67	1808.52	1371	567.361	1382.98	1995.77	1475.15	1715.16	1954	0	0	0	393.108	92.1369	0	1059.61	0	1587.63	0
	Cluster_100	0	1787.4	1288.79	1809.2	0	203.804	0.617404	243.918	1800.29	20.4127	1880.79	138.017	1865.87	1764.04	1910.84	1358.87	772.325	303.706	1850.54	1766.15	1926.4	1617.81	1813.5	1841.46	1040.42	2278.71	168.862	1397.44	2174	2286.63	1696.9	2039.48	1967.8	1835.68	660.908	2168.23	1472.08	3.34218	121.598	1951.97	2047.38	0	2.9218	0	0	1564.98	671.498	2041.98	1974.21	1165.22	1948.25	8.5865	0	1844.31	1016.6	0	1969.93	1360.33	989.876	0	83.3107	0	243.825	0	0	813.337	0	943.528	928.974	2022.98	1420.61	0	2065.53	0	2068.82	0	2116.92	1506.41	24.604	1278.64	1680.97	5.96794	2258.69	1655.07	1793.17	2281.66	536.34	1865.22	2099.05	1991.72	1887.63	1818.03	11.5384	444.321	1.74742	596.965	33.7661	2.56909	1334.86	0	2327.32	0
	Cluster_1000	0	1726.58	4340.43	1776.75	0	385.693	0.667669	326.695	1688.15	22.0745	1764.36	74.6268	1568.36	2170.37	1778.29	1330.68	1300.99	0	1547.67	1438.77	1794.79	1305.89	1846.12	2062.19	970.905	2278.08	0	1634.22	1978.07	2136.38	1895.21	1928.29	1958.73	1740.69	640.428	2005.58	1305.73	0	0	2023.48	1854.55	0	2.36976	21.4998	0	1870.29	726.167	2198.71	1758.71	0	1909.65	9.28556	4.60705	1681.68	0	63.4522	1877.75	1219.81	1213.19	0	0	0	0	10.2996	0	0	0	382.629	1268.97	1952.3	1442.94	0	2104.25	0	3355.87	0	1990.95	0	13.3035	1286.38	1773.84	0	1831.94	1282.26	1662.13	2202.62	473.921	1834.25	1680.4	1292.32	1747.84	1713.83	28.0751	240.247	0	542.947	10.9545	4.16737	1745.36	0	1510.08	0
	Cluster_1001	0	1742.53	1174.7	1638.48	0	293.861	0.333834	228.687	1700.23	0	1472.75	0	1685.3	36.3748	0	1543.99	626.402	0	1726.25	896.073	0	1476.94	1739.73	615.115	1277.16	2473.09	0	1476.07	1909.16	550.479	1624.47	0	29.0182	56.7885	670.925	788.959	0	0	0	1672.22	0	0	0	21.4998	0	1240.45	0	2002.4	1845.89	0	7.48884	4.64278	4.60705	1274.09	3298.09	126.904	828.279	0	749.325	12.1692	270.28	0	0	0	0	0	0	765.258	475.865	1744.61	1550.62	0	1379.41	17.2866	0	0	130.479	0	0	1179.93	1524.63	0	74.2677	1469.25	1569.79	389.456	1.72964	1730.15	0	1292.32	1667.52	1706.85	9.35837	840.865	1.41726	489.771	87.6362	6.25105	984.225	0	47.1899	0
	Cluster_1002	0	0	0	0	0	1597.87	0	3.62995	3.01994	0	0	0	0	0	0	0	0	0	0	0	0	4.17216	0	0	104.258	0	0	0	0	0	0	0	0	0	0	0	0	0	197.247	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
	Cluster_1003	1495.37	7.28202	80.0097	52.0906	0	1217.79	1224.14	1119.55	57.6446	1297.34	22.1566	1049.61	0	0	0	183.688	338.855	494.928	25.6292	126.792	0	209.574	57.7779	8.89156	510.612	0	963.142	52.9608	12.2166	55.9612	0	3.72184	0	0	904.993	0	139.264	1383.41	990.801	54.6732	0	0	1340.35	1619.95	1528.12	101.434	1094.29	62.7525	2.33556	0	77.116	1497.23	1272.8	55.9718	0	191.238	9.45583	55.0766	824.489	1662.67	1493.42	1295.29	2384.07	1660.73	1601.43	0	1089.06	768.801	557.746	0	151.453	1191.34	3.71528	1467.48	0	1445.92	0	0	1283.05	306.197	0	1397.23	49.741	13.4187	22.0877	7.77336	835.807	0	4.43092	0	14.7934	9.56537	1335.04	0	1460.84	922.221	1122.53	1549.07	197.756	1202.1	71.1126	1372.99
	...
	```
	
	* File name: `$OUTPUT_DIR/finalized/$BASENAME_proteinfamilies_nrm.tsv`
	* This file includes the normalized abundance of each protein family across samples (TSV format file).
	* Protein family abundance is reported in copies per million (CPM) units, which is "total sum scaling (TSS)"-style normalization: each sample is constrained to sum to 1 million. First, each protein family is normalized to RPK (reads per kilobase) units for gene length normalization; RPK units reflect relative gene (or transcript) copy number in the community. Then, RPK values are further sum-normalized (CPM) to adjust for differences in sequencing depth across samples. Further information can refer to the normalization approach in [HUMAnN](https://github.com/biobakery/humann). 
	
	**5. Clustering information for protein families**
	
	```
	>PRISM_7861_211956;Cluster_1;length=85;size=24;cluster=1
	PRISM_7861_211956
	PRISM_7791_142577
	PRISM_8794_132719
	...
	>PRISM_7855_54982;Cluster_2;length=559;size=12;cluster=2
	PRISM_7855_54982
	PRISM_8776_12844
	PRISM_8467_139927
	...
	```
	
	* File name: `$OUTPUT_DIR/finalized/$BASENAME_proteinfamilies.clstr`
	* This file includes the clustering information for protein families, formatted using an extension-fasta style based on the version of CD-hit clustering file (extension-fasta format file).
	
	**6. Sequences of protein families**  
	
	* File name: `$OUTPUT_DIR/finalized/$BASENAME_proteinfamilies.centroid.faa`
	* This file is the protein sequence for representatives of protein families (Fasta format file).

	**7. Intermediate output files**
	
	* Clustering results
		* MetaWIBELE clusters all representative sequences of gene families into protein families. 
		* All intermediate results are in the folder `$OUTPUT_DIR/clustering/`.
		
	* Global-homology based search results
		* MetaWIBELE queries each sequence in protein families against the UniRef90 database by performing protein-level search.
		* All intermediate results are in the folder `$OUTPUT_DIR/global_homology_annotation/`.
	
	* Domain-motif annotation results
		* MetaWIBELE uses the local-homology approach (domain/motif search) to characterize the secondary structures of protein families.
		* All intermediate results are in the folder `$OUTPUT_DIR/domain_motif_annotation/`.
	
	* Abundance-based annotation results
		* MetaWIBELE implements non-homology based strategy compromising (i) taxonomic annotation with phylogenetic binning, (ii) abundance profiling for protein families, and (iii) association with environmental parameters or phenotypes based on differential abundance. 
		* All intermediate results are in the folder `$OUTPUT_DIR/abundance_annotation/`.
	

#### MetaWIBELE-prioritize
* ##### Input files for prioritization
	* annotation file produced by MetaWIBELE-characterize workflow (TSV format file), e.g. [demo\_proteinfamilies_annotation.tsv](https://github.com/biobakery/metawibele/raw/master/examples/output/characterization/finalized/demo_proteinfamilies_annotation.tsv.gz)
	* annotation attribute file produced by MetaWIBELE-characterize workflow (TSV format file), e.g. [demo\_proteinfamilies\_annotation.attribute.tsv](https://github.com/biobakery/metawibele/raw/master/examples/output/characterization/finalized/demo_proteinfamilies_annotation.attribute.tsv.gz)

* ##### MetaWIBELE-prioritize workflow

	`$ metawibele prioritize --input-annotation $INPUT_ANNOTATION --input-attribute $INPUT_ATTRIBUTE --output $OUTPUT_DIR`
	
	* In the command replaces `$INPUT_ANNOTATION`, `$INPUT_ATTRIBUTE` with annotation and attribute files produced by MetaWIBELE-characterize workflow, `$OUTPUT_DIR` with the path to the folder to write output files. 
	* The workflow runs with the default settings for all main tool subtasks. These settings will work for most data sets. However, if you need to customize your workflow settings for the prioritization workflow to determine the optimum setting. Then apply these settings by using options for each task. You can customize your own configuration file.
		* For example, `--prioritization-config $myconfig_file` will modify the default settings when running the prioritization tasks.
	
* ##### Demo run of MetaWIBELE-prioritize

	`$ metawibele prioritize --input-annotation demo_proteinfamilies_annotation.tsv --input-attribute demo_proteinfamilies_annotation.attribute.tsv --output $OUTPUT_DIR`

* ##### Output files of MetaWIBELE-prioritize
	**1. unsupervised prioritization**
	
	```
	TID	familyID	evidence	value	rank	description	note
	1	Cluster_2	DNA_abundance	3976.775288	0.997474747	ranking based on single evidence	
	2	Cluster_2	DNA_prevalence	0.980392157	0.998737374	ranking based on single evidence	
	3	Cluster_2	priority_score	0.998105661	0.998105661	meta ranking based on multiple evidences	
	4	Cluster_1385	DNA_abundance	7579.266327	0.999368687	ranking based on single evidence	
	5	Cluster_1385	DNA_prevalence	0.960784314	0.994318182	ranking based on single evidence	
	6	Cluster_1385	priority_score	0.996837037	0.996837037	meta ranking based on multiple evidences	 
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
	TID	familyID	evidence	value	rank	description	note
	1	Cluster_459	DNA_within_phenotype_abundance	2006.049068	0.942838793	ranking based on single evidence	CD_vs_nonIBD
	2	Cluster_459	DNA_within_phenotype_prevalence	1	1	ranking based on single evidence	CD_vs_nonIBD
	3	Cluster_459	MaAsLin2_DA__mean_log	-3.55885575	0.892742453	ranking based on single evidence	CD_vs_nonIBD
	4	Cluster_459	MaAsLin2_DA__qvalue	1.00E-05	0.942122186	ranking based on single evidence	CD_vs_nonIBD
	5	Cluster_459	priority_score	0.942906085	0.942906085	meta ranking based on multiple evidences	CD_vs_nonIBD
	6	Cluster_1254	DNA_within_phenotype_abundance	1716.567921	0.846499679	ranking based on single evidence	CD_vs_nonIBD
	7	Cluster_1254	DNA_within_phenotype_prevalence	1	1	ranking based on single evidence	CD_vs_nonIBD
	8	Cluster_1254	MaAsLin2_DA__mean_log	-3.944110406	0.97430957	ranking based on single evidence	CD_vs_nonIBD
	9	Cluster_1254	MaAsLin2_DA__qvalue	7.56E-06	0.954662379	ranking based on single evidence	CD_vs_nonIBD
	10	Cluster_1254	priority_score	0.940027663	0.940027663	meta ranking based on multiple evidences	CD_vs_nonIBD
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
                     [--extension {.fastq.gz,.fastq}]
                     [--gene-call-type {prokka,prodigal,both}]
                     [--bypass-assembly] [--bypass-gene-calling]
                     [--bypass-gene-catalog]
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
                        provide the extension for paired fastq files using comma to separate, e.g. .R1.fastq.gz,.R2.fastq.gz | .R1.fastq,.R2.fastq
  --extension {.fastq.gz,.fastq}
                        provide the extension for all fastq files
                        [default: .fastq.gz]
  --gene-call-type {prokka,prodigal,both}
                        specify which type of gene calls will be used
                        [default: prodigal]
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
                        [default: /srv/export/hutlab11_nobackup/share_root/users/yancong/metawibele_demo/output/characterization/finalized]
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
                        [default: serial_requeue,serial_requeue,240]
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
* QC'ed shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format), e.g. "raw_reads" folder including:
	- [sample1_R1.fastq.gz](https://github.com/biobakery/metawibele/raw/master/examples/raw_reads/sample1_R1.fastq.gz)
	- [sample1_R2.fastq.gz](https://github.com/biobakery/metawibele/raw/master/examples/raw_reads/sample1_R2.fastq.gz)
	- [sample2_R1.fastq.gz](https://github.com/biobakery/metawibele/raw/master/examples/raw_reads/sample2_R1.fastq.gz)
	- [sample2_R2.fastq.gz](https://github.com/biobakery/metawibele/raw/master/examples/raw_reads/sample2_R2.fastq.gz)
* See the section on parallelization options to optimize the workflow run based on your computing resources. 
* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to modify the default settings, you can change the parameter settings.
	* For example, `--extension-paired "$R1_suffix,$R2_suffix"`, `--extension "$fastq_suffix"` (what are the following part after `$SAMPLE` in the input file names) will modify the default settings when running the assembly task.

#### Demo run of preprocessing workflow

`$ metawibele preprocess --input raw_reads/ --output $OUTPUT_DIR/ --output-basename demo --extension-paired "_R1.fastq.gz,_R2.fastq.gz" --extension ".fastq.gz"`

#### Output files of preprocessing workflow
**Main output files**

The following are the two main output files of the preprocessing utility that are used for MetaWIBELE:

```
$OUTPUT_DIR/finalized/$BASENMAE_genecatalogs.centroid.faa
$OUTPUT_DIR/finalized/$BASENMAE_genecatalogs_counts.all.tsv

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



**Other output files** 

**1. assembly results**
	
* `$OUTPUT_DIR/finalized/$BASENMAE_contig_sequence.fasta`: contig sequences (Fasta format file).
* The intermediate assembly outputs for each sample are in the `$OUTPUT_DIR/assembly/` folder.
	
**2. gene-calling results**
	
* `$OUTPUT_DIR/finalized/$BASENMAE_gene_info.tsv`: all gene calls information (TSV format file).
* `$OUTPUT_DIR/finalized/$BASENMAE_combined_gene_protein_coding.complete.sorted.fna`: nucleotide sequences for all complete ORFs sorted by gene length (Fasta format file).
* `$OUTPUT_DIR/finalized/$BASENMAE_combined_protein.complete.sorted.faa`: protein sequences for all complete ORFs sorted by protein length (Fasta format file).
* `$OUTPUT_DIR/finalized/$BASENMAE_combined_gene_protein_coding.sorted.fna`: nucleotide sequences for all ORFs (including partial genes) sorted by gene length (Fasta format file).
* `$OUTPUT_DIR/finalized/$BASENMAE_combined_protein.sorted.faa`: protein sequences for all ORFs (including partial genes) sorted by protein length.
* The intermediate gene-calling outputs from prodigal are in the `$OUTPUT_DIR/gene_calls/` folder. 
* The intermediate gene-annotation outputs from prokka are in the `$OUTPUT_DIR/gene_annotation/` folder.
	
**3. gene families**
	
* `$OUTPUT_DIR/finalized/$BASENMAE_genecatalogs.clstr`: clustering information for non-redundant gene families (extension-fasta format file).
* `$OUTPUT_DIR/finalized/$BASENMAE_genecatalogs_counts.all.tsv`: reads counts of gene families per sample (TSV format file).
* `$OUTPUT_DIR/finalized/$BASENMAE_genecatalogs.centroid.fna`: nucleotide sequences of representatives for gene families (Fasta format file).
* `$OUTPUT_DIR/finalized/$BASENMAE_genecatalogs.centroid.faa`: protein sequences of representatives for gene families.
* The intermediate mapping outputs for each sample are in the `$OUTPUT_DIR/mapping/` folder. 
	
----


# MetaWIBELE History #

## v0.4.6 2022-08-31 ##
* Added function to handle outputs of MSPminer >v1.0.0
* Printed the contents of the variable mysterr for debugging the dependency installs 
* The fixed version has been tagged v0.4.6

## v0.4.5 2022-06-27 ##
* Tweaked task names for the preprocessing workflow
* Fixed issue when running gridable jobs with one local job 
* Updated publication information
* The fixed version has been tagged v0.4.5

## v0.4.4 2021-11-24 ##
* Release v0.4.4
* Updated URLs of dependencies
* Added case when skipping mspminer-based taxonomic assignment

## v0.4.3 2021-11-22 ##
* Release v0.4.3
* Fixed maaslin2 incompatibility

## v0.4.2 2021-10-29 ##
* Release new version
* Fixed small issue on psortb annotator
* Cut down some space of demo data

## v0.4.1 2021-05-16 ##
* Fixed issues when running interproscan with all analyses
* Updated requirements for install
* Added a new module to check install
* Fixed issue when running utilities to prepare databases
* Fixed issue in config annotation types from interproscan  
* Keep temp directory
* Updated description texts about configs
* Added final taxonomy file with stratified taxon info
* Fixed small issue in utilities

## v0.4.0 2021-05-05 ##
* Added "--tshld-diff" "--tshld-classified" setting to global config file
* Added "--min_variance" "--reference" for MaAsLin2 >= v1.5.1
* Added "--global-config" option for taking the user provided global config file
* Added time stamp for each processing step
* Added default setting for uniref database and mutiple options for setting up the DB
* Added option for downloading demo database

## v0.3.9 2021-03-18 ##
* Added options for msp-based taxonomy assignment "--tshld-diff" "--tshld-classified"
* Renamed the name of config option for taxonomy assignment
* Fixed option setting for interproscan
* Added option for FDR correction
* Fixed issue when preprocessing contigs
* Fixed warnings when summarizing interproscan output

## v0.3.8 2020-10-02 ##
* Added bowtie2 options "--very-sensitive"
* Fixed issue of summarizing uniref annotations and contig sequences

## v0.3.7 2020-08-25 ##
* Added functions for downloading databases
* Added binary filtering for unsupervised prioritization
* Added config for interproscan executable file
* Fixed issue when formatting assembled protein sequences

## v0.3.6 2020-07-30 ##
* Make input/output options available for metawibele on the command line instead of the config file
* Add options for dependent files for each individual step
* Report all DDIs and microbe-human DDIs
* Fixed errors when running DDI-based annotations

## v0.3.5 2020-07-17 ##
* Skip phenotypic association by default
* Skip supervised prioritization if not associating with environmental/phenotypic parameters
* Turn off input/output options in workflows and only use info from the configuration file

## v0.3.4 2020-07-13 ##
* Update DDIs which are supported by PDB
* Update default settings of Maaslin2
* Update function to build dependent databases
* Set cases where running DDIs without domain annotations
* Add bypass modes for each individual domain/motif annotation step
* Rename finalized taxonomic annotation file

## v0.3.3 2020-06-18 ##
* Add utilities to build dependent uniref and domain databases
* Simplify the configurations in metawibele.cfg
* Deprecate sample_list and gene_catalog configs in metawibele.cfg
* Replace 'metawibele_workflow' with 'metawibele'
* Update configuration of depedent databases

## v0.3.2 2020-06-11 ##
* Use UniRef/UniProt databases in HUMAnN
* Remove to map assembled proteins to UniRef50
* Pack default dependent databases: domain databases, taxonomy databases
* Add function to download gloabl and local configuration files
* Remove specific utility scripts in the common folder
* Add option for readling compressed databases
* Balance gridable tasks
* Add option to preprocess outputs from only prokka, only prodigal, or both

## v0.3.1 2020-06-01 ##
* Convert input sources for gridable task
* Reformat annotations of PSORTb, InterProScan
* Add EffectiveT3 to predict effector (TBD)
* Make config file for MSPminer as optional

## v0.3 2020-05-21 ##
* Update different cases for binary filtering
* Update local and global configuration files
* Add Rscripts folder

## v0.2.2 2020-05-06 ##
* Update binary filtering approaches for supervised prioritization

## v0.2.1 2020-03-09 ##
* Add setup.py

## v0.2 2020-01-26 ##
* Add preprocessing workflow
* Update linear models for association analysis
* Modify abundance features

## v0.1 2019-12-16 ##
* Initial characterization and prioritization processing workflows added

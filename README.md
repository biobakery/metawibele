# MetaWIBELE User Manual

**MetaWIBELE** (**W**orkflow to **I**dentify novel **B**ioactive **Ele**ments in the microbiome) is a computational pipeline that identifies novel bioactive microbial gene products from metagenomes and finds new immunomodulatory gene families, especially targeting secreted/extracellular proteins to enrich for likely host interactors. The prioritized list of gene products can be further used for downstream experimental validation. MetaWIBELE is available as module of bioBakery [bioBakery repository](https://github.com/biobakery/metawibele).

## Citing MetaWIBELE

**A manuscript describing MetaWIBELE is currently in prep:**

Identifying Novel Bioactive Microbial Gene Products in Inflammatory Bowel Disease

**And feel free to link to MetaWIBELE in your Methods:**

[http://huttenhower.sph.harvard.edu/MetaWIBELE](http://huttenhower.sph.harvard.edu/MetaWIBELE)

**For additional information, read the** [MetaWIBELE Tutorial](https://github.com/biobakery/metawibele)

We provide support for MetaWIBELE users via our Google group. Please feel free to send any questions to the group by posting directly or emailing `<metawibele-users@googlegroups.com>`.
***

## Contents ##

* [Workflow](#workflow)
* [Install MetaWIBELE](#install-metawibele)
    * 1. [Requirements](#1.-requirements)
    * 2. [Installment](#2.-Installment)
    	* 2.1 [Download MetaWIBELE](#2.1-download-metawibele)
    	* 2.2 [Install MetaWIBELE](#2.2-install-metawibele)
    	* 2.3 [Download database](#2.3-download-databases)
    * 3. [Configuration for MetaWIBELE](#3.-configuration-for-metawibele)
    	* 3.1 [Global configration file](#3.1-global-configration-file)
    	* 3.2 [Configuration for characterization](#3.2-configuration-for-characterization)
    	* 3.3 [Configuration for prioritization](#3.3-configuration-for-prioritization)
    * 4. [Test with a demo run](#4.-test-with-a-demo-run) 
* [Quick-start Guide](#quick-start-guide)
    * 1. [How to run](#1.-how-to-run)
    * 2. [Standard Workflows](#2.-standard-workflows)
    	* 2.1 [MetaWIBELE-characterize workflow](#2.1-metawibele-characterize-workflow)
    		* [Input files for characterization](#input-files-for-characterization)
    		* [MetaWIBELE-characterize workflow](#metawibele-characterize-workflow)
    		* [To run a demo for characterization](#to-run-a-demo-for-characterization)
    		* [Output files for characterization](#output-files-for-characterization)
    	* 2.2 [MetaWIBELE-prioritize workflow](#2.2-metawibele-prioritize-workflow)
    		* [Input files for prioritization](#input-files-for-prioritization)
    		* [MetaWIBELE-prioritize workflow](#metawibele-prioritize-workflow)
    		* [To run a demo for prioritization](#to-run-a-demo-for-prioritization)
    		* [Output files for prioritization](#output-files-for-prioritization)
* [Guides to MetaWIBELE Utilities](#guides-to-metawibele-utilities)
	* [Qulity control for raw sequencing reads](#qulity-control-for-raw-sequencing-reads)
		* [Specific options for QC workflow](#specific-options-for-qc-workflow)
		* [How to run QC workflow](#how-to-run-qc-workflow)
		* [Example for running QC workflow](#example-for-running-qc-workflow) 
	* [Preprocessing sequencing reads into gene catalogs](#preprocessing-sequencing-reads-into-gene-catalogs)
		* [Specific options for preprocessing workflow](#specific-options-for-preprocessing-workflow)
		* [How to run preprocessing workflow](#how-to-run-preprocessing-workflow)
		* [Example for running preprocessing workflow](#example-for-running-preprocessing-workflow) 
		* [Output files of preprocessing workflow](#output-files-of-preprocessing-workflow)
* [Download MetaWIBELE resources](#download-metawibele-resources)
	* [Information of gene catalogs](#information-of-gene-catalogs)
	* [Characterization of protein families](#characterization-of-protein-families)
	* [Prioritization of protein families](#prioritization-of-protein-families)
    	
***


## Workflow
![workflow.png](https://www.dropbox.com/s/9wj1tufchuzgdmt/workflow.png?raw=1)
***


## Install MetaWIBELE
### 1. Requirements
```
* Python 3+ (tested with 3.7.4)
* Diamond (tested with v0.9.5)
* CD-hit (teset with version 4.7)
* Interproscan (tested with v5.31-70)
* Signalp (tested with v4.1)
* TMHMM (tested with v2.0c)
* Phobius (tested with 1.01)
* Psortb (tested with v3.0)
* MSPminer (tested with v2)
* MaAsLin2 (tested with version 0.99.12)
* AnADAMA2 (tested with version 0.5.0-devel)

(Optional)
* MEGAHIT (tested with v1.1.3)
* Prokka (tested with version 1.14-dev)
* Prodigal (tested with version 2.6)
* USEARCH (tested with version 9.0.2132_i86linux64)
* Bowtie2 (tested with version 2.3.2)
* SAMtools (tested with version 1.5)
* featureCounts (tested with version 1.6.2)

```

### 2. Installment
#### 2.1 Download MetaWIBELE
You can download the latest MetaWIBELE release or the development version. The source contains example files.

Option 1: Latest Release (Recommended)

* download [MetaWIBELE.zip](http://huttenhower.sph.harvard.edu/MetaWIBELE/MetaWIBELE.zip) and unpack the latested release of MetaWIBELE

Option 2: Development Version

* Create a clone of the Git repository 
	* `$ git clone https://github.com/biobakery/metawibele.git`
* You can always update to the latest version of the repository with:
	* `$ git pull --update`

#### 2.2 Install MetaWIBELE

* Move to the MetaWIBELE directory
	* `$ cd $MetaWIBELE_PATH`

* Install MetaWIBELE package
	* `$ python setup.py install`
	* If you do not have write permissions to '/usr/lib/', then add the option --user to the install command. This will install the python package into subdirectories of '~/.local'. Please note when using the "--user" install option on some platforms, you might need to add '~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `metawibele_workflow: command not found` when trying to run MetaWIBELE after installing with the '--user' option. Similarly, you can also specify the installment directory using '--prefix' option which will install the python package into the directory $YOUR\_INSTALL\_DIR specified using '--prefix'. You might need to add $YOUR\_INSTALL\_DIR/bin to your $PATH as it might not be included by default.

#### 2.3. Download databases
MetaWIBELE requires several databases that are needed to put in the **MetaWIBELE directory**. The versions used in the MetaWIBELE publication are available for download here. You need to download, unpack and put these databases in `$MetaWIBELE_PATH/data`.

* all databases: [metawibele_databases.tar](http://huttenhower.sph.harvard.edu/xxx) (58 GB)
	* UniRef databases: 
		* [uniref90.fasta.gz](http://huttenhower.sph.harvard.edu/xxx) (19 GB)
		* [uniref90.fasta.dmnd.gz](http://huttenhower.sph.harvard.edu/xxx)(20 GB)
		* [uniref50.fasta.gz](http://huttenhower.sph.harvard.edu/xxx) (6.4 GB)
		* [uniref50.fasta.dmnd.gz](http://huttenhower.sph.harvard.edu/xxx)(6.7 GB)
		* [uniref90.ann.tsv.gz](http://huttenhower.sph.harvard.edu/xxx)(5.2 GB)
		* [map\_UniRef90_UniRef50.dat.gz](http://huttenhower.sph.harvard.edu/xxx) (525 MB)
		* [BGC\_genes_unirefID.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (95 MB)
	* UniProt databases: 
		* [uniprot_taxonomy.map.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (43 MB)
		* [uniprot\_taxaID_mammalia.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (25 KB)
		* [uniprot\_taxaID_bac-arc-vir.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (1.6 MB)
		* [uniprot\_human_pfam.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (2.2 MB)
	* Pfam databases:
		* [Pfam_ann.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (214 KB)
		* [Pfam2GO.txt.gz](http://huttenhower.sph.harvard.edu/xxx) (98 KB)
	* DOMINE database:
		* [INTERACTION.txt.gz](http://huttenhower.sph.harvard.edu/xxx) (136 KB)
	* PDB database:
		* [pdb\_chain_taxonomy.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (27 MB)
		* [pdb\_chain_pfam.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (2.2 MB)
	* Expression Atlas databases:
		* [32\_Uhlen_Lab\_colon\_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (184 KB)
		* [Encode\_sigmoid_colon.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (128 KB)
		* [FANTOM5\_colon_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (147 KB)
		* [GTEx\_sigmoid_transverse\_colon.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (184 KB)
		* [Human\_protein_Atlas\_colon\_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (33 KB)
		* [Human\_proteome_map\_colon\_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (92 KB)
		* [Illumina\_Body_Map\_colon.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (133 KB)
	* example vignettes:
		* [vignettes_proteins.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (2.7 KB)

	
### 3. Configuration for MetaWIBELE
####3.1 Global configration file
When running MetaWIBELE, one global configuation file `metawibele.cfg` is required. You can copy `metawibele.cfg` file from the MetaWIBELE installment directory to your working directory, and then modify your configuration based on your datasets. You may need to specify:
* Input files:

```
[input]
study = my_study_name
metadata = /my/path/for/metadata/file/my_metadata.tsv
# if you have corresponding metadata file for metatranscriptomes (optional), otherwise specify as none
rna_metadata = none
sample_list = /my/path/for/sample/list/my_samples.tsv
gene_catalog = /my/path/for/gene_catalog/clusters/my_genecatalogs.clstr
gene_catalog_nuc = /my/path/for/gene_catalog/nucleotide/sequences/my_genecatalogs.centroid.fna
gene_catalog_prot = /my/path/for/gene_catalog/protein/sequences/	my_genecatalogs.centroid.faa
gene_catalog_count = /my/path/for/gene_catalog/reads/counts/my_genecatalogs.counts.all.tsv
# if you have corresponding reads counts table for metatranscriptomes (optional), otherwise specify as none
gene_catalog_rna_count = none
``` 

* Output files:
	
```
[output]
basename = my_output_prefix
working_dir = /my/workding/direcoty/path/
```

* parameter setting for MaAsLin2:

```
[maaslin2]
maaslin2_output = my_maaslin2_output_directory
maaslin2_cmmd = /my/path/Maaslin2/R/Maaslin2.R
fixed_effects = metadata2,metadata3,metadata4
random_effects = metadata5
# specify nested effect if you have, e.g. disease activity within diagnosis (diagnosis:CD,UC), otherwise specify it as none
nested_effects = metadata1:category1,category2
maaslin2_cores = 25
tshld_prevalence = 0.10
tshld_qvalue = 0.05
effect_size = coef
abundance_detection_level = 0 
min_abundance = -20 
min_prevalence = 0 
max_significance = 0.05
normalization = NONE
transform = NONE
analysis_method = LM
plot_heatmap = FALSE
plot_scatter = FALSE
# specify contrast status for the metadata type: commas used for multiple values for the same metadata, semicolons used for multuple metadatas
contrast_status = metadata2:contrastCategory1,contrastCategory2
# specify reference status for the metadata type: commas used for multiple values for the same metadata, semicolons used for multuple metadatas
ref_status = metadata2:contrastCategory1_vs_refCategory1,contrastCategory2_vs_refCategory2
```

####3.2 Configuration for characterization
You can specify the characterization modules in this configuration file. Default will run all modules for characterization. 
	
`my_characterization.cfg`:
	
```
[protein_family]
#protein family annotation based on global similarity: [yes] process this step, [no] skip this step
uniref = yes
#BGC annotation based on global similarity: [yes] process this step, [no] skip this step
antismash = yes

[domain_motif]
#domain annotation: [yes] process this step, [no] skip this step
interproscan = yes
#Pfam2GO to annotate GOs: [yes] process this step, [no] skip this step
pfam2go = yes
#domain-domain interaction from DOMINE database: [yes] process this step, [no] skip this step
domine = yes
#DDI with SIFTS evidence: [yes] process this step, [no] skip this step
sifts = yes
#DDI with human expression from ExpAtlas database: [yes] process this step, [no] skip this step
expatlas = yes
#subcellular annotation: [yes] process this step, [no] skip this step
psortb = yes

[abundance]
#summary DNA abundance: [label] provide label for DNA abundance, e.g. DNA_abundance, [no] skip this step
dna_abundance = DNA_abundance
#summary RNA abundance: [label] provide label for RNA abundance, e.g. RNA_abundance, [no] skip this step
rna_abundance = RNA_abundance
#differential abundance based on DNA abundance: [label] provide label for DA annotation, e.g. MaAsLin2_DA, [no] skip this step
dna_da = MaAsLin2_DA
#differential abundance based on RNA abundance: [label] provide label for DE annotation, e.g. MaAsLin2_DE, [no] skip this step
rna_de = MaAsLin2_DE-RNA

[integration]
#summarize annotation info: [yes] process this step, [no] skip this step
summary_ann = yes
#generate finalized annotations: [yes] process this step, [no] skip this step
finalization = yes
``` 

####3.3 Configuration for prioritization
You can specify the prioritization modules in this config file. Default will run all configurations of prioritization in MetaWIBELE installed directory. 

`my_prioritization.cfg`:
	
```
[unsupervised]
#threshold for highly prioritized list based on top % prioritized list: [proportion] top percentage, NaN means ignoring
tshld_priority = 0.25 
#threshold for highly prioritized list based on priority score: [score] score threshold, NaN means ignoring
tshld_priority_score = NaN
#weight parameter of prevalence: [proportion] values of parameter
beta = 0.50
#weight value of prevalence for calculating priority: [proportion] weight value
DNA-all-nonIBD_prevalence = 0.50
#weight value of mean abundance for calculating priority: [proportion] weight value
DNA-all-nonIBD_abundance = 0.50

[supervised]
#threshold for highly prioritized list based on top % prioritized list: [proportion] top percentage, NaN means ignoring
tshld_priority = 0.50 
#threshold for highly prioritized list based on priority score: [score] score threshold, NaN means ignoring
tshld_priority_score = NaN
#ecological property (mean abundance) for calculating priority: [required] required item, [optional] optional item, [none] ignoring
#MaAsLin2_DA__mean_prevalent_abundance = required
DNA-all-within-phenotype_abundance = required
#ecological property (prevalence) for calculating priority: [required] required item, [optional] optional item, [none] ignoring
DNA-all-within-phenotype_prevalence = required
#association property (q values) for calculating priority: [required] required item, [optional] optional item, [none] ignoring
MaAsLin2_DA__qvalue = required
#association property (fold change) for calculating priority: [required] required item, [optional] optional item, [none] ignoring
MaAsLin2_DA__foldChange = required

##Binary filtering for selection subset
#All [vignette_type], [cluster_file] items should be true: [vignette_type] required interested function type; [cluster_file] required subset of proteins
#All [required] items should be true: [required] required item 
#At least one [optional] item should be true: [optional] optional item
#All [none] items will be ignored: [none] ignoring
[filtering]
#interested functional vignettes type: [vignette_type] vignettes types, e.g. pilin | superfamily | ect.
#vignettes = molybdopterin
# interested protein list file: [cluster_file] proteins file, specific_cluster.tsv
#clusters = specific_cluster.tsv
#biochemical annotation for filtering: [required] required item, [optional] optional item, [none] ignoring
ExpAtlas_interaction = required
#association annotation for filtering: [required] required item, [optional] optional item, [none] ignoring
#MaAsLin2_DE-RNA = required
#population property (novel taxonomy) for filtering: [required] required item, [optional] optional item, [none] ignoring
#Novel_taxonomy = required
#population property (accessory gene) for filtering: [required] required item, [optional] optional item, [none] ignoring 
#MSPminer_accessory = required
#property of (multiple MSPs) for filtering: [required] required item, [optional] optional item, [none] ignoring
#MSPminer_multi-msp = required
DOMINE_interaction = none
SIFTS_interaction = none
Denovo_signaling = optional
PSORTb_extracellular = optional
PSORTb_cellWall = optional
PSORTb_outerMembrane = optional
Denovo_transmembrane = optional
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

##Application: categorial modules
[moduling]
#categorical clusters: [required] required item, [optional] optional item, [none] ignoring 
#MSPminer_family-module = required
MSPminer_module = required
#specific clusters: [cluster_file] specific proteins file, specific_cluster.tsv
#clusters = structural_clusters.tsv
```

### 4. Test with a demo run
* Run with default configuration
	* First, run MetaWIBELE for characterization:
`$ metawibele_workflow characterize --input examples/input/ --output examples/`
	
	* Then, run MetaWIBELE for prioritization based on the output of the above characterization step:
`$ metawibele_workflow prioritize --input examples/input/ --output examples/`

* Run with specified configuration (see more details about configuration files for each workflow in the [Configuration for MetaWIBELE](#markdown-header-3.-configuration-for-metawibele) session)
	* Run MetaWIBELE for characterization:
`$ metawibele_workflow characterize --characterization-config my_characterization.cfg --input examples/input/ --output examples/`
	
	* Run MetaWIBELE for prioritization based on the output of the above characterization step:
`$ metawibele_workflow prioritize --prioritization-config my_prioritization.cfg --input examples/input/ --output examples/`

***


## Quick-start Guide
### 1. How to run
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


### 2. Standard Workflows
#### 2.1 MetaWIBELE-characterize workflow
* #####Input files for for characterization
	* clustering information for non-redudant gene catalogs using extended-fasta format, e.g. [demo_genecatalogs.clstr]()
	* nucleotide sequences for non-redudant gene catalogs, e.g. [demo_genecatalogs.centroid.fna]()
	* protein seqeuences for non-redudant gene catalogs, e.g. [demo_genecatalogs.centroid.faa]()
	* reads counts table for non-redudant gene catalogs, e.g. [demo\_genecatalogs_counts.all.tsv]()
	* metadata file, e.g. [demo\_mgx_metadata.tsv]()
	* sample list, e.g. [demo\_MGX_samples.tsv]()
	* all the above information can be specified in the `metawibele.cfg` file.

* #####MetaWIBELE-characterize workflow
	`$ metawibele_workflow characterize --input $INPUT --output $OUTPUT`
	* In the command replace $INPUT with the path to the folder containing your fastq input files and $OUTPUT with the path to the folder to write output files. See the section on parallelization options to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings to run all modules. These settings will work for most data sets. However, if you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting. You can specify which modules you want to run in your own configuration file.
	* For example, `--characterization-config=$myconfig_file` will modify the default settings when running the characterization modules.
	
* #####To run a demo for characterization

	`$ metawibele_workflow characterize --characterization-config my_characterization.cfg --input examples/input/ --output examples/`

* #####Output files for characterization
	**1. Annotation file**
	
	```
	familyID    annotation  feature category    method  AID
Cluster_1   demo    study   project Shotgun NA
Cluster_1   Cluster_1   protein_family  Denovo_clustering   CD-hit  Cluster_1__Denovo_clustering
Cluster_1   UniRef90_A0A3E2UKI3 strong_homology UniRef90_homology   UniRef90    Cluster_1__UniRef90_homology
Cluster_1   Faecalibacterium prausnitzii    Species Taxonomy_characterization   Taxonomy_annotation Cluster_1__Taxonomy_characterization
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
14  Cluster_1__Taxonomy_characterization    taxa_id 853
15  Cluster_1__Taxonomy_characterization    taxa_lineage    k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii
16  Cluster_1__Taxonomy_characterization    LCA_Tax Faecalibacterium prausnitzii
17  Cluster_1__Taxonomy_characterization    LCA_TaxID   853
18  Cluster_1__Taxonomy_characterization    Rep_Tax Faecalibacterium prausnitzii
19  Cluster_1__Taxonomy_characterization    Rep_TaxID   853
20  Cluster_1__Taxonomy_characterization    msp_name    msp_unknown
21  Cluster_1__Taxonomy_characterization    msp_taxa_name   NA
22  Cluster_1__Taxonomy_characterization    msp_taxa_id NA
23  Cluster_1__DNA-CD_abundance mean_abundance  145.88748347267776
24  Cluster_1__DNA-CD_abundance mean_prevalent_abundance    164.03938233794182
25  Cluster_1__DNA-CD_abundance prevalence  0.889344262295082
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
	* Protein family abundance is reported in copies per million (CPM) units, which is "total sum scaling (TSS)"-style normalization: each sample is constrained to sum to 1 million. First, each protein family is normalized to RPK (reads per kilobase) units for gene length normalization; RPK units reflect relative gene (or transcript) copy number in the community. Then, RPK values are further sum-normalized (CPM) to adjust for differences in sequencing depth across samples. Further information can refer to the normalization approach in [HUMAnN2] (). 
	
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
	

#### 2.2 MetaWIBELE-prioritize workflow
* #####Input files for prioritiztion
	* anntation file produced by MetaWIBELE-characterize workflow, e.g. [demo_proteinfamilies_annotation.tsv]()
	* attribute file produced by MetaWIBELE-characterize workflow, e.g. [demo_proteinfamilies_annotation.attribute.tsv]()
	* all the above information can be specified in the `prioritization.cfg` file.

* #####MetaWIBELE-prioritize workflow

	`$ metawibele_workflow prioritize --input $INPUT --output $OUTPUT`
	* In the command replace $INPUT with the path to the folder containing your fastq input files and $OUTPUT with the path to the folder to write output files. See the section on parallelization options to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings for all main tool subtasks. These settings will work for most data sets. However, if you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting. Then apply these settings by using options for each task. You can specify your own configuration file.
	* For example, `--prioritization-config=$myconfig_file` will modify the default settings when running the prioritization tasks.
	
* #####To run a demo for prioritization

	`$ metawibele_workflow prioritize --prioritization-config my_prioritization.cfg --input examples/characterization/ --output examples/`

* #####Output files for prioritization
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
	$OUTPUT\_DIR/prioritization/$BASENAME\_unsupervised\_prioritization.priority.tsv
	* These files are the results of unsupervised prioritization based on ecological properties. Each of protein family has a numeric priority score.
	* `$BASENAME_unsupervised_prioritization.rank.tsv` is the overall ranking for all protein families.
	* `$BASENAME_unsupervised_prioritization.priority.tsv`is one subset of unsupervised prioritization results with higher priority scores.
	
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
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.priority.tsv
	* These file are the results of supervised prioritization by combing ecological properties and assoiciation with host phenotypes. Each of protein family has a numeric priority score.
	* `$BASENAME_supervised_prioritization.rank.tsv` is the overall ranking for all protein families.
	* `$BASENAME_supervised_prioritization.priority.tsv`is one subset of supervised prioritization results with higher priority scores.

	**3. supervised prioritization: binary filtering**
	
	```
	familyID    DNA_within_phenotype_abundance__value   DNA_within_phenotype_abundance__percentile  DNA_within_phenotype_prevalence__value  DNA_within_phenotype_prevalence__percentile MaAsLin2_DA__coef__value    MaAsLin2_DA__coef__percentile   MaAsLin2_DA__qvalue__value  MaAsLin2_DA__qvalue__percentile priority_score
Cluster_196|CD.dysbiosis_vs_CD.non_dysbiosis    214.2212814279175   0.9686878727634195  0.9816933638443935  0.9924410075341273  -105.831926383983   0.9866302186878728  8.01016319636424e-09    0.991574918607252   0.9847393815443188
Cluster_25293|CD.dysbiosis_vs_CD.non_dysbiosis  120.48241777803203  0.933051689860835   0.9931350114416476  0.9973394335728671  -67.2740317585006   0.9606858846918489  7.98001508992818e-11    0.9997266197778164  0.9719079557712381
Cluster_129|CD.dysbiosis_vs_CD.non_dysbiosis    344.4962334050344   0.9944333996023856  0.9725400457665904  0.9875923117089788  -76.7386486741326   0.9705765407554672  9.37199044695546e-06    0.9324005268782464  0.9706435542254798
Cluster_101|CD.dysbiosis_vs_CD.non_dysbiosis    535.9237779633867   0.9983598409542743  0.9633867276887872  0.9799338588159237  -96.7983809913901   0.9830019880715706  2.91093846771617e-05    0.9144816959514874  0.9678369794322969
Cluster_273|CD.dysbiosis_vs_CD.non_dysbiosis    99.57769844622432   0.9164512922465209  0.9694835680751174  0.9847576895342766  -62.7194765911449   0.9544234592445328  1.06787781457237e-08    0.9904316922235753  0.960601551965445
Cluster_533|CD.dysbiosis_vs_CD.non_dysbiosis    222.18690086498862  0.9715208747514911  0.9084507042253521  0.9206803093219286  -65.7308323652607   0.9584493041749503  0.000216885848009021    0.8523995327683476  0.9233657031704845
...
	```
	* File name: $OUTPUT_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.filter.tsv
	* This file is the results of supervised filtering of protein families based on binary annotation features.
	* The default filtering approach is requiring that each of prioritized protein family should: 
		*  be annotated to domain-domain interaction with host
		*  have at least one of the following features: signaling, extracellular, cellWall, outerMembrane, transmembrane 
	
	**4. finalized prioritization**
	
	```
	TID familyID evidence    value   rank    description note
1   Cluster_14393   DNA_within_phenotype_abundance  844.0252184037556   0.9995526838966203	ranking based on single evidence  CD.dysbiosis_vs_CD.non_dysbiosis
2   Cluster_14393   DNA_within_phenotype_prevalence 0.971830985915493   0.986622572543949	ranking based on single evidence   CD.dysbiosis_vs_CD.non_dysbiosis
3   Cluster_14393   MaAsLin2_DA__coef   -480.492828810768   0.9995526838966203	ranking based on single evidence  CD.dysbiosis_vs_CD.non_dysbiosis
4   Cluster_14393   MaAsLin2_DA__qvalue 4.31866766640033e-11    1.0	ranking based on single evidence CD.dysbiosis_vs_CD.non_dysbiosis
5   Cluster_14393   priority_score  0.9963995495816218  0.9963995495816218  meta ranking based on multiple evidences    CD.dysbiosis_vs_CD.non_dysbiosis
	```
	* File name: 
		$OUTPUT\_DIR/prioritization/$BASENAME\_unsupervised\_prioritization.rank.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_unsupervised\_prioritization.priority.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.priority.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.filter.table.tsv
	* These file are formated the prioritization results in the same way.

***


## Guides to MetaWIBELE Utilities
###Qulity control for raw sequencing reads
A utitlity workflow in MetaWIBELE package for QC, which integrates KneadData quality control pipeline (http://huttenhower.sph.harvard.edu/KneadData) with additional automatic adapter detection, trimming low-quality read bases and removing potentially conmanitationed reads.

#### Specific options for QC workflow
```
usage: metawibele_qc_workflow [-h] [--version] --sample-file SAMPLE_FILE
                      [--trimmomatic-options TRIMMOMATIC_OPTIONS]
                      [--additional-options ADDITIONAL_OPTIONS]
                      [--remove-intermediate-output REMOVE_INTERMEDIATE_OUTPUT]
                      [--contaminant-db CONTAMINANT_DB]
                      [--file-extension FILE_EXTENSION] [--threads THREADS] -o
                      OUTPUT [-i INPUT] [--local-jobs JOBS]
                      [--grid-jobs GRID_JOBS] [--grid GRID]
                      [--grid-partition GRID_PARTITION]
                      [--grid-benchmark {on,off}]
                      [--grid-options GRID_OPTIONS]
                      [--grid-environment GRID_ENVIRONMENT] [--dry-run]
                      [--skip-nothing] [--quit-early]
                      [--until-task UNTIL_TASK] [--exclude-task EXCLUDE_TASK]
                      [--target TARGET] [--exclude-target EXCLUDE_TARGET]
                      [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

A workflow to run kneaddata on the input files provided to perform quality control.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --sample-file SAMPLE_FILE
                        Sample files including sample names (string).
  --trimmomatic-options TRIMMOMATIC_OPTIONS
                        options for trimmomatic (string): ILLUMINACLIP:/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 (optional)
                        [default: none]
  --additional-options ADDITIONAL_OPTIONS
                        additional_options (string): Additional options when running kneaddata (optional)
                        [default: none]
  --remove-intermediate-output REMOVE_INTERMEDIATE_OUTPUT
                        remove_intermediate_output (bool): Remove intermediate output files.
                        [default: True]
  --contaminant-db CONTAMINANT_DB
                        Select reference sequences for the contamination you are trying to remove. It is KneadData databases including the indexed redernece sequences.
                        [default: none]
  --file-extension FILE_EXTENSION
                        Extension of input fastq files (string)
                        [default: .R1.fastq.gz,.R2.fastq.gz]
  --threads THREADS     number of threads/cores for each task to use
                        [default: 6]
  -o OUTPUT, --output OUTPUT
                        Write output to this directory
  -i INPUT, --input INPUT
                        Find inputs in this directory 
                        [default: /n/home00/yancong/projects/R24_HMBR/src/assembly-based/metawibele/tools]
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
                        [default: serial_requeue,general,240]
  --grid-benchmark {on,off}
                        Benchmark gridable tasks 
                        [default: on]
  --grid-options GRID_OPTIONS
                        Grid specific options that will be applied to each grid task
  --grid-environment GRID_ENVIRONMENT
                        Commands that will be run before each grid task to set up environment
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

* `--sample-file`: a txt file listing sample names or identifiers corresponding to sequences, e.g. [sample.txt]()
* `--input`: the input directory where a set of fastq (or fastq.gz) files (single-end or paired-end) are stored. The files are expected to be named `$SAMPLE.fastq.gz`,`$SAMPLE.R1.fastq.gz`, or `$SAMPLE.R2.fastq.gz` where `$SAMPLE` is the sample name or identifier corresponding to the sequences. `$SAMPLE` can contain any characters except spaces or periods.
* `--contaminant-db`: select reference sequences for the contamination you are trying to remove. It is KneadData databases including the indexed redernece sequences. See more details in [KneadData muanual](https://github.com/biobakery/kneaddata).
* `--file-extension`: the extension for fastq (or fastq.gz) files. It will be specified as **".R1.fastq.gz,R2.fastq.gz"** if the files named as `$SAMPLE.R1.fastq.gz` and `$SAMPLE.R2.fastq.gz`.
* `--output`: the ouput directory. QC'ed reads will generated into the subfolder `$OUTPUT_DIR/kneaddata/`

#### How to run QC workflow
* See the section on parallelization options to optimize the workflow run based on your computing resources. 
* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting, you can change the parameter settings.
* For example, `--file-extension="$R1_suffix,$R2_suffix"` (what are the follwong part after `$SAMPLE` in the input file names) will modify the default settings when running the assembly task.

#### Example for running QC workflow

`$ metawibele_qc_workflow --thread 10 --sample-file sample.txt --input examples/paired/ --output examples/cleaned_reads/ --grid-jobs 10 --grid-partition serial_requeue`

---

###Preprocessing sequencing reads into gene catalogs
A utility workflow in MetaWIBELE package for preprocessing metagenomes reads, used for (i) metagenomic assembly, (ii) open reading frame prediction, (iii) non-redundant gene catalogs construction and (iv) gene abundance estimation.

#### Specific options for preprocessing workflow
```
usage: metawibele_preprocessing_workflow [-h] [--version] [--threads THREADS]
                                 --sample-file SAMPLE_FILE --extension-paired
                                 EXTENSION_PAIRED
                                 [--extension-orphan EXTENSION_ORPHAN]
                                 [--assembly ASSEMBLY]
                                 [--gene-calling GENE_CALLING]
                                 [--gene-catalog GENE_CATALOG] -o OUTPUT
                                 [-i INPUT] [--local-jobs JOBS]
                                 [--grid-jobs GRID_JOBS] [--grid GRID]
                                 [--grid-partition GRID_PARTITION]
                                 [--grid-benchmark {on,off}]
                                 [--grid-options GRID_OPTIONS]
                                 [--grid-environment GRID_ENVIRONMENT]
                                 [--dry-run] [--skip-nothing] [--quit-early]
                                 [--until-task UNTIL_TASK]
                                 [--exclude-task EXCLUDE_TASK]
                                 [--target TARGET]
                                 [--exclude-target EXCLUDE_TARGET]
                                 [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

A workflow to preprocess shotgun sequencing reads of metagenomes with tasks of metagenomic assembly, gene calling, building gene catalogs and generating gene abundance for each sample.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --threads THREADS     number of threads/cores for each task to use
                        [default: 20]
  --sample-file SAMPLE_FILE
                         the list file of samples
  --extension-paired EXTENSION_PAIRED
                        indicates the extension for paired fastq files, e.g. R1.fastq.gz,R2.fastq.gz
  --extension-orphan EXTENSION_ORPHAN
                        indicates the extension for orphan fastq files
                        [default: none]
  --assembly ASSEMBLY   indicates whether or not do assembly
                        [default: True]
  --gene-calling GENE_CALLING
                        indicates whether or not call ORFs
                        [default: True]
  --gene-catalog GENE_CATALOG
                        indicates whether or not build gene catalogs
                        [default: True]
  -o OUTPUT, --output OUTPUT
                        Write output to this directory
  -i INPUT, --input INPUT
                        Find inputs in this directory 
                        [default: /n/home00/yancong/projects/R24_HMBR/src/assembly-based/metawibele/tools]
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
                        [default: serial_requeue,general,240]
  --grid-benchmark {on,off}
                        Benchmark gridable tasks 
                        [default: on]
  --grid-options GRID_OPTIONS
                        Grid specific options that will be applied to each grid task
  --grid-environment GRID_ENVIRONMENT
                        Commands that will be run before each grid task to set up environment
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

* `--sample-file`: a txt file listing sample names or identifiers corresponding to sequences, e.g. [sample.txt]()
* `--input`: the input directory where a set of fastq (or fastq.gz) files (single-end or paired-end) passing through QC are stored. The files are expected to be named `$SAMPLE.paired_R1.gz`, `$SAMPLE.paired_R2.gz`, `$SAMPLE.orphan_R1.gz` and `$SAMPLE.orphan_R2.gz` where `$SAMPLE` is the sample name or identifier corresponding to the sequences. `$SAMPLE` can contain any characters except spaces or periods.
* `--extension-paired` indicates the extension for paired fastq files. It should be specified as **".paired_R1.fastq.gz,.paired_R2.fastq.gz"** if the paired fastq files are `$SAMPLE.paired_R1.gz` and `$SAMPLE.paired_R2.gz`  
* `--extension-orphan` indicates the extension for orphan fastq files. It should be specified as **".orphan_R1.fastq.gz,.orphan_R2.fastq.gz"** if the orphan fastq files are `$SAMPLE.orphan_R1.gz` and `$SAMPLE.orphan_R2.gz`  
* `--output`: the ouput directory. 

#### How to run preprocessing workflow
* See the section on parallelization options to optimize the workflow run based on your computing resources. 
* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting, you can change the parameter settings.
* For example, `--extension-paried="$R1_suffix,$R2_suffix"`, `--extension-orphan="$orphan_suffix"` (what are the follwong part after `$SAMPLE` in the input file names) will modify the default settings when running the assembly task.

#### Example for running preprocessing workflow

`$ preprocessing_workflow --thread 10 --sample-file sample.txt --input examples/paired/ --output examples/preprocessing/ --extension-paired ".paired_R1.fastq.gz,.paired_R2.fastq.gz" --extension-orphan ".orphan_R1.fastq.gz,.orphan_R2.fastq.gz" --grid-jobs 10 --grid-partition serial_requeue`

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
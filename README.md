# MetaWIBELE User Manual

**MetaWIBELE** (**W**orkflow to **I**dentify novel **B**ioactive **Ele**ments in the microbiome) is a computational pipeline that identifies novel bioactive microbial gene products from metagenomes and finds new immunomodulatory gene families, especially targeting secreted/extracellular proteins to enrich for likely host interactors. The prioritized list of gene products can be further used for downstream experimental validation. MetaWIBELE is available as module of bioBakery [bioBakery repository](https://github.com/biobakery).

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
    * [Requirements](#requirements)
    * [Installation](#installation)
    	* [Download MetaWIBELE](#download-metawibele)
    	* [Install MetaWIBELE](#install-metawibele)
    	* [Download database](#download-databases)
    * [Configuration for MetaWIBELE](#configuration-for-metawibele)
    	* [Global configration file](#global-configration-file)
    	* [Configuration for characterization](#configuration-for-characterization)
    	* [Configuration for prioritization](#configuration-for-prioritization)
    * [Test with a demo run](#test-with-a-demo-run) 
* [Quick-start Guide](#quick-start-guide)
    * [How to run](#how-to-run)
    * [Standard Workflows](#standard-workflows)
    	* [MetaWIBELE-characterize workflow](#metawibele-characterize-workflow)
    		* [Input files for characterization](#input-files-for-characterization)
    		* [MetaWIBELE-characterize workflow](#metawibele-characterize-workflow)
    		* [To run a demo for characterization](#to-run-a-demo-for-characterization)
    		* [Output files for characterization](#output-files-for-characterization)
    	* [MetaWIBELE-prioritize workflow](#metawibele-prioritize-workflow)
    		* [Input files for prioritization](#input-files-for-prioritization)
    		* [MetaWIBELE-prioritize workflow](#metawibele-prioritize-workflow)
    		* [To run a demo for prioritization](#to-run-a-demo-for-prioritization)
    		* [Output files for prioritization](#output-files-for-prioritization)
* [Guides to MetaWIBELE Utilities](#guides-to-metawibele-utilities)
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
### Requirements
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

### Installation
#### Download MetaWIBELE
You can download the latest MetaWIBELE release or the development version. The source contains example files.

Option 1: Latest Release (Recommended)

* download [MetaWIBELE.zip](http://huttenhower.sph.harvard.edu/MetaWIBELE/MetaWIBELE.zip) and unpack the latested release of MetaWIBELE

Option 2: Development Version

* Create a clone of the Git repository 
	* `$ git clone https://github.com/biobakery/metawibele.git`
* You can always update to the latest version of the repository with:
	* `$ git pull --update`

#### Install MetaWIBELE

* Move to the MetaWIBELE directory
	* `$ cd $MetaWIBELE_PATH`

* Install MetaWIBELE package
	* `$ python setup.py install`
	* If you do not have write permissions to '/usr/lib/', then add the option --user to the install command. This will install the python package into subdirectories of '\~/.local/'. Please note when using the '--user' install option on some platforms, you might need to add '\~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `metawibele_workflow: command not found` when trying to run MetaWIBELE after installing with the '--user' option. 
	* Similarly, you can also specify the installation directory using '--prefix' option which will install the python package into the directory $YOUR\_INSTALL\_DIR that is a directory on PYTHONPATH and which Python reads ".pth" files from. You might need to add $YOUR\_INSTALL\_DIR/bin to your $PATH as it might not be included by default.

#### Download databases
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

	
### Configuration for MetaWIBELE
#### Global configration file
When running MetaWIBELE, one global configuation file `metawibele.cfg` is required. You can copy `metawibele.cfg` file from the MetaWIBELE installation directory to your working directory, and then modify your configuration based on your datasets. You may need to specify:

* Input files:

```
[input]
study = my_study_name
metadata = /my/path/for/metadata/file/my_metadata.tsv
sample_list = /my/path/for/sample/list/my_samples.tsv
gene_catalog = /my/path/for/gene_catalog/clusters/my_genecatalogs.clstr
gene_catalog_prot = /my/path/for/gene_catalog/protein/sequences/	my_genecatalogs.centroid.faa
gene_catalog_count = /my/path/for/gene_catalog/reads/counts/my_genecatalogs.counts.all.tsv

``` 

* Output files:
	
```
[output]
basename = my_output_prefix
working_dir = /my/workding/direcoty/path/
```

* Computational resources

```
[computation]
# the number of cores that you’re requesting
threads = 10
# the amount of memory (in MB) that you will be using for your job 
memory = 20000
```

* Abundance configs

```
[abundance]
# binning gene catalogs based on co-abundance info and do taxonomic annotation: [config] provide mspminer config file, [no] skip this step
mspminer = configs/MSPminer_setting.cfg 
# the method for normalization: [cpm]:copies per million units (sum to 1 million); [relab]: relative abundance (sum to 1) 
normalize = cpm
```

* Parameter settings for MaAsLin2:

```
[maaslin2]
maaslin2_output = my_maaslin2_output_directory
maaslin2_cmmd = /my/path/Maaslin2/R/Maaslin2.R
min_abundance = 0 
min_prevalence = 0 
max_significance = 0.05
normalization = NONE
transform = LOG
analysis_method = LM
fixed_effects = metadata1,metadata2,metadata3,metadata4
random_effects = metadata5
correction = BH
standardize = TRUE
plot_heatmap = FALSE
heatmap_first_n = FALSE
plot_scatter = FALSE
maaslin2_cores = 10
tshld_prevalence = 0.10
tshld_qvalue = 0.05
effect_size = mean(log)
abundance_detection_level = 0
# specify one main phenotype used for prioritization
phenotype = metadata2
# flag reference status as control for phenotype varibles; use semicolon to seperate variables
flag_ref = "metadata1:control_status1;metadata2:control_status2"
# specify case and control status paries for phenotype variables; use comma to seperate each comparison within the sample variable; use semicolon to seperate variables
case_control_status = "metadata1:case_status1|control_status1; metadata2:case_status2|control_status2,case_status3|control_status2"

```

#### Configuration for characterization
You can specify the characterization modules in this configuration file. Default will run all modules for characterization. 
	
`my_characterization.cfg`:
	
```
[protein_family]
#protein family annotation based on global similarity: [yes] process this step, [no] skip this step
uniref = yes

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
#differential abundance based on DNA abundance: [label] provide label for DA annotation, e.g. MaAsLin2_DA, [no] skip this step
dna_da = MaAsLin2_DA

[integration]
#summarize annotation info: [yes] process this step, [no] skip this step
summary_ann = yes
#generate finalized annotations: [yes] process this step, [no] skip this step
finalization = yes
``` 

#### Configuration for prioritization
You can specify the prioritization modules in this config file. Default will run all configurations of prioritization in MetaWIBELE installed directory. 

`my_prioritization.cfg`:
	
```
## Mandatory ranking
[unsupervised]
# weight parameter of prevalence: [proportion] values of parameter
beta = 0.50
# weight value of prevalence for calculating priority: [proportion] weight value
DNA-nonIBD.non-dysbiosis_prevalence = 0.50
# weight value of mean abundance for calculating priority: [proportion] weight value
DNA-nonIBD.non-dysbiosis_abundance = 0.50

[supervised]
# ecological property (mean abundance) for calculating priority: [required] required item, [optional] optional item, [none] ignoring
#MaAsLin2_DA__mean_prevalent_abundance = required
DNA-within-phenotype_abundance = required
# ecological property (prevalence) for calculating priority: [required] required item, [optional] optional item, [none] ignoring
#MaAsLin2_DA__prevalence = required
DNA-within-phenotype_prevalence = required
# association property (q values) for calculating priority: [required] required item, [optional] optional item, [none] ignoring
MaAsLin2_DA__qvalue = required
# association property (coefficient) for calculating priority: [required] required item, [optional] optional item e.g. coef | mean(log) | log(FC), [none] ignoring
MaAsLin2_DA__mean(log) = required


##Binary filtering for selection subset
#All [vignette_type], [cluster_file] items should be true: [vignette_type] required interested function type; [cluster_file] required subset of proteins
#All [required] items should be true: [required] required item 
#At least one [optional] item should be true: [optional] optional item
#All [none] items will be ignored: [none] ignoring
[filtering]
# interested functional vignettes type: [vignette_type] vignettes types, e.g. pilin | superfamily | ect.
#vignettes = pilin

# significant associations for filtering: [required] required item, [optional] optional item, [none] ignoring
#MaAsLin2_DA-sig = required

# biochemical annotation for filtering: [required] required item, [optional] optional item, [none] ignoring
ExpAtlas_interaction = required
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

### Test with a demo run
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
	* clustering information for non-redudant gene catalogs using extended-fasta format, e.g. [demo_genecatalogs.clstr]()
	* nucleotide sequences for non-redudant gene catalogs, e.g. [demo_genecatalogs.centroid.fna]()
	* protein seqeuences for non-redudant gene catalogs, e.g. [demo_genecatalogs.centroid.faa]()
	* reads counts table for non-redudant gene catalogs, e.g. [demo\_genecatalogs_counts.all.tsv]()
	* metadata file, e.g. [demo\_mgx_metadata.tsv]()
	* sample list, e.g. [demo\_MGX_samples.tsv]()
	* all the above information can be specified in the `metawibele.cfg` file.

* ##### MetaWIBELE-characterize workflow
	`$ metawibele_workflow characterize --input $INPUT --output $OUTPUT`
	* In the command replace $INPUT with the path to the folder containing your fastq input files and $OUTPUT with the path to the folder to write output files. See the section on parallelization options to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings to run all modules. These settings will work for most data sets. However, if you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting. You can specify which modules you want to run in your own configuration file.
	* For example, `--characterization-config=$myconfig_file` will modify the default settings when running the characterization modules.
	
* ##### To run a demo for characterization

	`$ metawibele_workflow characterize --characterization-config my_characterization.cfg --input examples/input/ --output examples/`

* ##### Output files for characterization
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
	

#### MetaWIBELE-prioritize workflow
* ##### Input files for prioritiztion
	* anntation file produced by MetaWIBELE-characterize workflow, e.g. [demo_proteinfamilies_annotation.tsv]()
	* attribute file produced by MetaWIBELE-characterize workflow, e.g. [demo_proteinfamilies_annotation.attribute.tsv]()
	* all the above information can be specified in the `prioritization.cfg` file.

* ##### MetaWIBELE-prioritize workflow

	`$ metawibele_workflow prioritize --input $INPUT --output $OUTPUT`
	* In the command replace $INPUT with the path to the folder containing your fastq input files and $OUTPUT with the path to the folder to write output files. See the section on parallelization options to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings for all main tool subtasks. These settings will work for most data sets. However, if you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting. Then apply these settings by using options for each task. You can specify your own configuration file.
	* For example, `--prioritization-config=$myconfig_file` will modify the default settings when running the prioritization tasks.
	
* ##### To run a demo for prioritization

	`$ metawibele_workflow prioritize --prioritization-config my_prioritization.cfg --input examples/characterization/ --output examples/`

* ##### Output files for prioritization
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
	
	***3.1 Select interested subset based on specific functions***
	
	Setting `my_prioritization.cfg` as following:
	
	```
	##Binary filtering for selection subset
	# All [vignette_type], [cluster_file] items should be true: [vignette_type] required interested function type; [cluster_file] required subset of proteins
	# All [required] items should be true: [required] required item 
	# At least one [optional] item should be true: [optional] optional item
	# All [none] items will be ignored: [none] ignoring
	[filtering]
	# interested functional vignettes type: [vignette_type] vignettes types, e.g. pilin | superfamily | ect.
	vignettes = pilin

	# significant associations for filtering: [required] required item, [optional] optional item, [none] ignoring
	#MaAsLin2_DA-sig = required

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
	
	* Use provied the file containing the interested functions which is formated as this example file: [vignettes_proteins.tsv](). By defualt, MetaWIBELE will use the file provided by the package. Users can also specify their own file for filtering specific functions with the required formmat by using `--interested-function` parameter when run MetaWIBELE prioritization workflow. E.g. `$ metawibele_workflow prioritize --prioritization-config my_prioritization.cfg --interested-function my_own_interested_functions_file --input examples/characterization/ --output examples/` 
	* File name: $OUTPUT_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.selected.tsv
	* This file is the results of supervised filtering of protein families based on pilin related functions.
	
	***3.2 Select interested subset annotated with at least one of specific biochemical annotations***
	 
	Setting `my_prioritization.cfg` as following:
	
	```
	##Binary filtering for selection subset
	# All [vignette_type], [cluster_file] items should be true: [vignette_type] required interested function type; [cluster_file] required subset of proteins
	# All [required] items should be true: [required] required item 
	# At least one [optional] item should be true: [optional] optional item
	# All [none] items will be ignored: [none] ignoring
	[filtering]
	# interested functional vignettes type: [vignette_type] vignettes types, e.g. pilin | superfamily | ect.
	#vignettes = pilin

	# significant associations for filtering: [required] required item, [optional] optional item, [none] ignoring
	#MaAsLin2_DA-sig = required

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
	
	* File name: $OUTPUT_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.selected.tsv
	* This file is the results of supervised filtering of protein families based on biochemical annotations.
	* These settings require that each of prioritized protein family should 1) be annotated to domain-domain interaction with host, and 2) have at least one of the following features: signaling, extracellular, cellWall, outerMembrane, transmembrane 
	
    ***3.3. Select interested subset annotated with multiple specific biochemical annotations simultaneously***
	 
	 Setting `my_prioritization.cfg` as following:
	
	```
	##Binary filtering for selection subset
	# All [vignette_type], [cluster_file] items should be true: [vignette_type] required interested function type; [cluster_file] required subset of proteins
	# All [required] items should be true: [required] required item 
	# At least one [optional] item should be true: [optional] optional item
	# All [none] items will be ignored: [none] ignoring
	[filtering]
	# interested functional vignettes type: [vignette_type] vignettes types, e.g. pilin | superfamily | ect.
	#vignettes = pilin

	# significant associations for filtering: [required] required item, [optional] optional item, [none] ignoring
	#MaAsLin2_DA-sig = required

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
	
	 * File name: $OUTPUT_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.selected.tsv
	 * This file is the results of supervised filtering of protein families based on biochemical annotations.
	 * These settings require that each of prioritized protein family should 1) be annotated to domain-domain interaction with host and 3) predicted as signal peptides, and 3) have at least one of the following features: extracellular, cellWall, outerMembrane, transmembrane 
	
	
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
		$OUTPUT\_DIR/prioritization/$BASENAME\_unsupervised\_prioritization.priority.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.priority.table.tsv
		$OUTPUT\_DIR/prioritization/$BASENAME\_supervised\_prioritization.rank.filter.table.tsv
	* These file are formated the prioritization results in the same way.

***


## Guides to MetaWIBELE Utilities

### Preprocessing sequencing reads into gene catalogs
A utility workflow in MetaWIBELE package for preprocessing metagenomes reads, used for (i) metagenomic assembly, (ii) open reading frame prediction, (iii) non-redundant gene catalogs construction and (iv) gene abundance estimation.

#### Specific options for preprocessing workflow
```
usage: metawibele_workflow preprocess [-h] [--version] [--threads THREADS]
                     [--extension-paired EXTENSION_PAIRED]
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
                        [default: 20]
  --extension-paired EXTENSION_PAIRED
                        provide the extension for paired fastq files using comma to seperate, e.g. .R1.fastq.gz,.R2.fastq.gz | .R1.fastq,.R2.fastq
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
                        [default: metawibele]
  -o OUTPUT, --output OUTPUT
                        Write output to this directory
  -i INPUT, --input INPUT
                        Find inputs in this directory 
                        [default: /n/home00/yancong/lib/python/lib/python3.7/site-packages/metawibele-0.3-py3.7.egg/metawibele]
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
* `--extension-paired` indicates the extension for paired fastq files using comma to seperate. It should be specified as **".R1.fastq.gz,.R2.fastq.gz"** if the paired fastq files are `$SAMPLE.R1.fastq.gz` and `$SAMPLE.R2.fastq.gz`  
* `--extension` indicates the extension for all fastq files. It should be specified as **".fastq.gz"** if the fastq files are `$SAMPLE.fastq.gz` 
* `--output`: the ouput directory. 

#### How to run preprocessing workflow
* See the section on parallelization options to optimize the workflow run based on your computing resources. 
* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to determine the optimum seeting, you can change the parameter settings.
* For example, `--extension-paried="$R1_suffix,$R2_suffix"`, `--extension="$fastq_suffix"` (what are the follwong part after `$SAMPLE` in the input file names) will modify the default settings when running the assembly task.

#### Example for running preprocessing workflow

`$ preprocessing_workflow --input /my/path/preprocessing/cleaned_reads/ --output /my/path/preprocessing/ --extension-paired ".R1.fastq.gz,.R2.fastq.gz" --extension ".fastq" --output-basename demo --local-jobs 10`

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

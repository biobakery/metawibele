# MetaWIBELE User Manual

**MetaWIBELE** (**W**orkflow to **I**dentify novel **B**ioactive **Ele**ments in the microbiome) is a computational pipeline that identifies novel bioactive microbial gene products from metagenomes and finds new immunomodulatory gene families, especially targeting secreted/extracellular proteins to enrich for likely host interactors. The prioritized list of gene products can be further used for downstream experimental validation.


## Citing MetaWIBELE

**A manuscript describing MetaWIBELE is currently in prep:**

Identifying Novel Bioactive Microbial Gene Products in Inflammatory Bowel Disease

**And feel free to link to MetaWIBELE in your Methods:**

[http://huttenhower.sph.harvard.edu/MetaWIBELE](http://huttenhower.sph.harvard.edu/MetaWIBELE)

We provide support for MetaWIBELE users via our Google group. Please feel free to send any questions to the group by posting directly or emailing `<metawibele-users@googlegroups.com>`.
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

### 2. Install MetaWIBELE
Install the MetaWIBELE software and Python dependencies with `pip`:

```
$ pip install metawibele
```

You can also clone the MetaWIBELE package from bitbucket with git (`git`):

```
$ git clone xxx
```

Or download and install MetaWIBELE package directly:

```
1. download source codes
$ wget MetaWIBELE.zip
$ unzip MetaWIBELE.zip
$ cd $MetaWIBELE_PATH # Move to the MetaWIBELE directory

2. Install MetaWIBELE
$ python setup.py install 
# If you do not have write permissions to '/usr/lib/', then add the option --user to the install command. 
# This will install the python package into subdirectories of '~/.local'. 
# Please note when using the "--user" install option on some platforms, you might need to add '~/.local/bin/' to your $PATH as it might not be included by default. 
# You will know if it needs to be added if you see the following message metawibele.py: command not found when trying to run MetaWIBELE after installing with the '--user' option.

```

### 3. Download the databases
MetaWIBELE requires several databases that are needed to put in the `data/` folder in the MetaWIBELE directory. These databases are included in the MetaWIBELE source package. The versions used in the MetaWIBELE publication are available for download here:

* all databases: [data.tar.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
* UniRef databases: 
	* [uniref90.fasta.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [uniref90.fasta.dmnd.gz](http://huttenhower.sph.harvard.edu/xxx)(xx)
	* [uniref50.fasta.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [uniref50.fasta.dmnd.gz](http://huttenhower.sph.harvard.edu/xxx)(xx)
	* [uniref90.ann.tsv.gz](http://huttenhower.sph.harvard.edu/xxx)(xx)
	* [uniref50.ann.tsv.gz](http://huttenhower.sph.harvard.edu/xxx)(xx)
	* [map\_UniRef90_UniRef50.dat.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [BGC\_genes_unirefID.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
* UniProt databases: 
	* [uniprot_taxonomy.map.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx GB)
	* [uniprot\_taxaID_mammalia.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [uniprot\_taxaID_bac-arc-vir.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [uniprot\_human_pfam.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
* Pfam databases:
	* [Pfam_ann.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [Pfam2GO.txt.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
* DOMINE database:
	* [pdb\_chain_taxonomy.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
* Expression Atlas databases:
	* [32\_Uhlen_Lab\_colon\_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [Encode\_sigmoid_colon.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [FANTOM5\_colon_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [GTEx\_sigmoid_transverse\_colon.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [Human\_protein_Atlas\_colon\_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [Human\_proteome_map\_colon\_rectum.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	* [Illumina\_Body_Map\_colon.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
* Essential genes:
	* [DEG\_genes_unirefID.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
* BGC genes:
	* [BGC\_genes_unirefID.tsv.gz](http://huttenhower.sph.harvard.edu/xxx) (xx)
	
### 4. Test with a demo run

```
$ metawibele_workflow $WORKFLOW --input examples/ --output test_demo/
```
***


## Quick-start Guide
### 1. How to run
* Run **characterization** workflow

```
metawibele_workflow characterize \
	--config $characterization_conf \
	--input $INPUT \
	--output $OUTPUT
```

* Run **prioritization** workflow

```
metawibele_workflow prioritize \
	--quantitative-config $quantitative_conf \
	--selection-config $selective_conf \
 	--input $INPUT \
 	--output $OUTPUT 
```

### 2. Standard Workflows
![workflow.png](https://bitbucket.org/repo/ReAr48a/images/3226399479-workflow.png)


#### 2.1 MetaWIBELE-characterize workflow
* Inputs
* Configuration file
* To run a demo
* Outputs

#### 2.2 MetaWIBELE-prioritize workflow
* Inputs
* Configuration file
* To run a demo
* Outputs
***


## Guides to MetaWIBELE Utilities
[TBD]
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

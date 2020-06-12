#!/usr/bin/env Rscript
##########################################################################
# Function: Utilitis for MAasLin2
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 05/29/2019
##########################################################################

# get the options
args <- commandArgs(T)
action = args[1]
data_in = args[2]
data_out = args[3]
pcl_utils = args[4]
split_num = args[5]

source(pcl_utils)


##########################################
# Prepare reference status 
prepare_ref <- function(metafile, meta, outfile) {
  # extract metadata
  metadata <- read.table(metafile, row.names = 1, head=TRUE, sep="\t", quote="", comment.char = "", check.names = FALSE)
  metadata$diagnosis <- as.factor(metadata$diagnosis) # convert to category
  metadata$diagnosis <- relevel(metadata$diagnosis, ref = meta) # set 'nonIBD' as reference
  
  # build formula
  #formula <- as.formula(~ age + diagnosis + antibiotic + immunosuppressant + mesalamine + steroids)
  formula <- as.formula(paste("~ ", paste(colnames(metadata), collapse= "+")))
  
  # model.matrix
  mymatrix <- model.matrix(formula, metadata)
  mymatrix <- mymatrix[,-1]
  new_col_name <- c("ID", colnames(mymatrix)) 

  # output matrix
  write.table(t(new_col_name), file = outfile, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE) 
  write.table(mymatrix, file = outfile, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = TRUE, col.names = FALSE, quote=FALSE)
}


##########################################
# Prepare nested/interacted effect 
prepare_effect <- function(metafile, outfile) {
  # extract metadata
  metadata <- read.table(metafile, head=TRUE, sep="\t", quote="", comment.char = "", check.names = FALSE)
  attach(metadata)
  new <- interaction(metadata$diagnosis, metadata$active)
  metadata <- data.frame(metadata, diagnosis2=new)
  new_col_name <- c(colnames(metadata))
  
  # output metadata
  write.table(t(new_col_name), file = outfile, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE) 
  write.table(metadata, file = outfile, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
}


##########################################
# convert table
convert_table <- function(infile, outfile) {
  a <- read.table(infile)
  b <- t(a)
  write.table(b, file = outfile, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
}


##########################################
# split table
split_table <- function(infile, split_num, outfile) {
	a <- read.table(infile, row.names = 1, head=TRUE, sep="\t", quote="", comment.char = "", check.names = FALSE)
	rownum <- nrow(a)
	split_size <- as.integer(rownum/as.numeric(split_num))
	i <- 1
	j <- split_size
	if (j < i) {
		j <- rownum
	}
	mynum <- 1
	while (i <= rownum)
	{
		myfile <- paste(outfile, mynum, ".pcl", sep="")
		b <- a[i:j,]
  		new_col_name <- c("ID", colnames(a))
  		write.table(t(new_col_name), file = myfile, append = FALSE, sep = "\t", eol = "\n", na = "nan", row.names = FALSE, col.names = FALSE, quote=FALSE) 
  		#write.table(t(new_col_name), file = myfile, append = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE, quote=FALSE) 
  		write.table(b, file = myfile, append = TRUE, sep = "\t", eol = "\n", na = "nan", row.names = TRUE, col.names = FALSE, quote=FALSE)
		i <- j + 1
		j <- min(i + split_size - 1, rownum)
		mynum <- mynum + 1
	}
}


##########################################
# Filter out feature data
filter_feature <- function(pclfile, outfile) {
  pcl <- pcl.read(pclfile, metadata.rows = NA, rowfeatures=T, sep="\t")
  sds <- pcl.apply.f(pcl, sd(x, na.rm=T))
  pcl <- pcl %>% pcl.filter.f(keep=sds > median(sds) / 2)
  pcl.write(pcl, outfile, rowfeatures=T) 
}


##########################################
# Combine results and multiplicity correction
multiplicity_correction <- function (infile, outfile) {
	print(infile)
	paras <- read.table(infile, header=TRUE, sep="\t", quote="", comment.char = "", check.names = FALSE)
	new_col_name <- c(colnames(paras), "qvalue")

	# multiplicity correction based on all metadata
	paras$qvalue <- as.numeric(p.adjust(as.numeric(paras$pval), "BH"))
	write.table(t(new_col_name), file = outfile, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
	write.table(paras, file = outfile, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
	
	# multiplicity correction based on each variable
	outfile_per <- gsub(".tsv", ".correct_per_variable.tsv", outfile)
	paras_per <- ddply(paras, .(metadata), transform, qvalue=as.numeric(p.adjust(as.numeric(pval), "BH")))
	write.table(t(new_col_name), file = outfile_per, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
	write.table(paras_per, file = outfile_per, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
	
	# multiplicity correction based on each level
	outfile_per <- gsub(".tsv", ".correct_per_level.tsv", outfile)
	paras_per <- ddply(paras, .(metadata, value), transform, qvalue=as.numeric(p.adjust(as.numeric(pval), "BH")))
	write.table(t(new_col_name), file = outfile_per, append = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
	write.table(paras_per, file = outfile_per, append = TRUE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE, quote=FALSE)
}


###########################################
##################  main  #################
# filter out features with little variation
if (action == "filter")
{
	filter_feature (data_in, data_out)
}

# Multiplicity correction
if (action == "correct")
{
	multiplicity_correction (data_in, data_out)
}

# Modify metadata: add nested effect
if (action == "modification")
{
	prepare_effect (data_in, data_out)
}

# split pcl table
if (action == "split")
{
	split_table (data_in, split_num, data_out)
}

library(cluster)
library(gplots)
library(RColorBrewer)
library("gplots")
library("devtools")
#require("TROM")
#require("GMD")
#library("heatmap3")
#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
source("~/usr/heatmap.3.R")

#mycol <- colorRampPalette(brewer.pal(8, "RdBu"))(80)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
mycol <- rf(80)
#mycol <- c("#ffe874", "#d80000")

args <- commandArgs(T)
sample_file <- args[1]
cond_file <- args[2]
myname <- args[3]
sample_prefix <- args[4]


###### Condtion info ####
#mycolor <- c(CD="#8B1A1A", nonIBD="#56B4E9", UC="#E69F00")
#dysbiosis_colors <- c("CD.dysbiosis"="#8B1A1A", "CD.non_dysbiosis"="#b97575", nonIBD="#56B4E9", "UC.dysbiosis"="#E69F00", "UC.non_dysbiosis"="#d8ba77")
conds <- read.table(cond_file, header=T, sep="\t")
col.cell <- c("#8B1A1A", "#b97575", "#56B4E9", "#E69F00", "#d8ba77")[conds$Condition]
#col.cell = as.matrix(t(col.cell))
#colnames(col.cell) = c("Diagnosis")[conds$Condition]
col.cell <- as.matrix(col.cell)
rownames(col.cell) <- conds$Group


######## Clutser samples ########
### Prepare data ###
print('Reading matrix file.')
print(sample_file)
primary_data <- read.table(sample_file, header=T, com='', sep="\t", check.names=F)
mydata <- primary_data
rownames <- mydata$familyID
mydata <- mydata[, -1]
#mydata = mydata[, -1]
data <- as.matrix(mydata)
rownames(data) <- rownames


# convert data 
data <- log2(data + 0.000001)
gene_cor <- NULL
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
myheatcol <- mycol
data <- t(scale(t(data), scale=F)) # normalized rows, mean substracted

# do clustering based on euclidean distance
dist_no_na <- function(mat) {
	edist <- dist(mat, method='euclidean')
	edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1
	edist[is.na(edist)] <- 0
	edist[is.nan(edist)] <- 0
	return(edist)
}

# redo the sample clustering according to the normalized distance values.
sample_dist <- dist(t(data), method='euclidean')
sub1 <- data[, substr(colnames(data), 1, 12) == "CD.dysbiosis"]
sub2 <- data[, substr(colnames(data), 1, 16) == "CD.non_dysbiosis"]
sub3 <- data[, substr(colnames(data), 1, 12) == "UC.dysbiosis"]
sub4 <- data[, substr(colnames(data), 1, 16) == "UC.non_dysbiosis"]
sub5 <- data[, substr(colnames(data), 1, 6)  == "nonIBD"]
#sample_dist1 = dist(t(sub1), method='euclidean')
sample_dist1 <- dist_no_na (t(sub1))
sample_dist2 <- dist_no_na (t(sub2))
sample_dist3 <- dist_no_na (t(sub3))
sample_dist4 <- dist_no_na (t(sub4))
sample_dist5 <- dist_no_na (t(sub5))

# redo the sample clustering with complete method. Complete linkage method for hierarchical clustering by default. This particular clustering method defines the cluster distance between two clusters to be the maximum distance between their individual components
hc_samples <- hclust(sample_dist, method='complete')
hc_samples <- as.dendrogram(hc_samples)
#hc_samples = hclust(sample_dist, method='average')
hc_samples1 <- hclust(sample_dist1, method='complete')
hc_samples1 <- as.dendrogram(hc_samples1)
hc_samples2 <- hclust(sample_dist2, method='complete')
hc_samples2 <- as.dendrogram(hc_samples2)
hc_samples3 <- hclust(sample_dist3, method='complete')
hc_samples3 <- as.dendrogram(hc_samples3)
hc_samples4 <- hclust(sample_dist4, method='complete')
hc_samples4 <- as.dendrogram(hc_samples4)
hc_samples5 <- hclust(sample_dist5, method='complete')
hc_samples5 <- as.dendrogram(hc_samples5)
neworder1 <- colnames(sub1)[order.dendrogram(hc_samples1)]
neworder2 <- colnames(sub2)[order.dendrogram(hc_samples2)]
neworder3 <- colnames(sub3)[order.dendrogram(hc_samples3)]
neworder4 <- colnames(sub4)[order.dendrogram(hc_samples4)]
neworder5 <- colnames(sub5)[order.dendrogram(hc_samples5)]
#cn <- c(neworder1, neworder2, neworder3, neworder4, neworder5)
#cn <- colnames(data)[order.dendrogram(hc_samples)]
#heatmap_data <- data[,cn]
#col.cell <- as.matrix(col.cell[cn,])
heatmap_data <- data

# redo the gene clustering according to normalized distance values.
#gene_dist = dist(data, method='euclidean')
gene_dist = dist_no_na (heatmap_data)
# redo the gene clustering
hc_genes = hclust(gene_dist, method='complete')

mypdf = paste(sample_prefix, "pdf", sep = ".")
pdf(mypdf, width=28, height=22)
data_out = paste("Genes vs. samples in ", myname, sep="")
nc = length(colnames(heatmap_data))
nr = length(rownames(heatmap_data))

# dendrogram is computed and reordered based on row means and column means. In heatmaps.3, there is method parameter. The agglomeration method to be used by hclust function. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
heatmap.3(heatmap_data, col=myheatcol, dendrogram='both', scale="none", Rowv=as.dendrogram(hc_genes), Colv=hc_samples, ColSideColors=col.cell, ColSideColorsSize=1.5, RowSideColorsSize=1.2, density.info="none", trace="none", key=TRUE, cexRow = 1.8 + 1/log10(nr), cexCol = 1.8 + 1/log10(nc), KeyValueName="relative abundance", keysize=40, labCol=FALSE, margins=c(5, 30), lhei=c(1,6), lwid=c(1,3), main=data_out, cex.main=40)	# using hclust to cluster genes with method "average"
#heatmap.3(heatmap_data, col=myheatcol, dendrogram='both', scale="none", Rowv=as.dendrogram(hc_genes), Colv=hc_samples, ColSideColorsSize=1.5, RowSideColorsSize=1.2, density.info="none", trace="none", key=TRUE, cexRow = 1.8 + 1/log10(nr), cexCol = 1.8 + 1/log10(nc), KeyValueName="relative abundance", keysize=20, labCol=FALSE, margins=c(0, 20), lhei=c(1,6), lwid=c(1,3), main=data_out, cex.main=40)	# using hclust to cluster genes with method "average"
#heatmap.3(heatmap_data, RowIndividualColors=myheatcol, dendrogram='row', scale="none", 
#          Rowv=as.dendrogram(hc_genes), Colv=NULL, ColSideColors=col.cell, ColSideColorsSize=1.5,
#          RowSideColorsSize=1.2, density.info="none", trace="none", key=TRUE, 
#          cexRow = 1.8 + 1/log10(nr), cexCol = 1.8 + 1/log10(nc), 
#          KeyValueName="relative abundance", keysize=20, labCol=FALSE, margins=c(0, 20), 
#          lhei=c(1,6), lwid=c(1,3), main=data_out, cex.main=40) 
legend("topright", legend=c("CD.dysbiosis", "CD.non_dysbiosis", "UC.dysbiosis", "UC.non_dysbiosis", "nonIBD"), 
       fill=c("#8B1A1A", "#b97575", "#E69F00", "#f0c566", "#56B4E9"), 
       border=FALSE, bty="n", y.intersp=0.7, cex=2)
dev.off()

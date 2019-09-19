#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
#  PCL I/O

library(plyr)
library(dplyr)

pcl.read <- function(datafile, metadata.rows = NA, rowfeatures=T, sep="\t") {
    # Read a pcl structure from a file
    # if metadata.rows is NA, it tries to guess the number of rows automatically
    # if metadata.rows is " ", a blank line is expected after the metadata

    library(data.table)
    starttime <- proc.time()["elapsed"]

    # Read in everything
    if (endsWith(datafile, ".gz")) {
        file <- gzfile(datafile, "r")
        raw <- read.delim(file, stringsAsFactors=F, fill=T, blank.lines.skip=F, sep=sep)
        close(file)
    } else {
        raw <- fread(datafile, data.table=F, stringsAsFactors=F, fill=T, sep=sep, header=T)
    }
    featnames <- raw[,1]
    dups <- duplicated(featnames)
    if (any(dups)) {
        featnames[dups] <- sprintf("%s_#%d", featnames[dups], cumsum(dups)[dups])
    }
    rownames(raw) <- featnames
    raw <- raw[,-1]

    if (rowfeatures) {
        # Transpose if rows are features in the datafile
        raw <- as.data.frame(t(raw), stringsAsFactors=F)
    }

    # Which columns are numeric?
    sampnames <- rownames(raw)
    if (any(colwise(mode)(raw) != "numeric")) {
        raw <- colwise(type.convert)(raw)
    }
    rownames(raw) <- sampnames

    if (is.na(metadata.rows)) {
        # Try to figure out how many metadata there are
        first_allnum <- suppressWarnings(max(which(!colwise(is.numeric)(raw)))) + 1
        first_k__ <- suppressWarnings(min(which(grepl("^k__", colnames(raw)))))
        first_pwy <- suppressWarnings(min(which(grepl("PWY", colnames(raw)))))
        first_nonchar <- min(which(grepl("[^a-zA-Z0-9_]", colnames(raw))))
        first_feat <- min(first_k__, first_pwy, first_nonchar)

        metadata.rows <- max(first_allnum, ifelse(is.finite(first_feat), first_feat, 1) - 1)
    } else if (is.character(metadata.rows) && metadata.rows == " ") {
        blank <- min(which(featnames==""))
        if (!is.finite(blank)) {
            stop(sprintf("metadata.rows is \" \", but no blank line was found in %s to denote the end of the metadata block.", datafile))
        }
        metadata.rows <- blank - 1
        featnames <- featnames[-blank]
        raw <- raw[,-blank]
    }

    # Pull out data and metadata
    metadata <- raw[,seq_len(metadata.rows)]
    x <- as.matrix(if (metadata.rows > 0) {
        raw[,-seq_len(metadata.rows)]
    } else {
        raw
    })

    # Report time taken
    runtime <- proc.time()["elapsed"] - starttime
    cat(sprintf("Loaded %s in %gs\n", datafile, runtime))

    return (list(x=x, meta=metadata, ns=dim(x)[1], nf=dim(x)[2]))
}

pcl.write <- function(dat, datafile, rowfeatures=T) {
    # Write a pcl structure to a file

    if (rowfeatures) {
        # Write the sample names
        write.table(cbind("ID", matrix(rownames(dat$x), nrow=1)), file=datafile, sep="\t",
                    quote=F, col.names=F, row.names=F)

        # Write the metadata
        if (dim(dat$meta)[2] > 0) {
            write.table(t(as.matrix(dat$meta)), file=datafile, sep="\t",
                        quote=F, append=T, col.names=F)
        }

        # Append the data
        write.table(t(dat$x), file=datafile, sep="\t",
                    quote=F, append=T, col.names=F)
    } else {
        # Merge into one big table
        write.table(dat$x, file=datafile, sep="\t",
                    quote=F, append=F, col.names=NA)
    }
}

pcl.empty <- list(nf=0, ns=0, x=matrix(0,0,0), meta=data.frame())

pcl.make <- function(x, meta, metaf) {
    x <- as.matrix(x)
    dat <- list(x=x, nf=dim(x)[2], ns=dim(x)[1])
    dat$meta <- if (missing(meta)) {
        data.frame(row.names=rownames(x))
    } else {
        stopifnot(nrow(meta) == nrow(x)) # Metadata does not have the same number of samples as the data
        meta[match(rownames(x), rownames(meta)),, drop=F]
    }
    if (!missing(metaf)) {
        stopifnot(nrow(metaf) == ncol(x)) # Feature metadata does not have the same number of features as the data
        dat$metaf <- metaf[match(colnames(x), rownames(metaf)),, drop=F]
    }
    return (dat)
}

# -----------------------------------------------------------------------------
#  Function application over samples and features

pcl.apply.s <- function(dat, expr, fn = substitute(expr), enclos = parent.frame()) {
    # Function application over samples

    # Todo: add metaf support

    evalfor <- function(i) {
        e <- c(as.list(dat$meta[i,,drop=F]), as.list(dat$x[i,,drop=F]))
        e$Name <- rownames(dat$x)[i]
        e$Index <- i
        e$x <- dat$x[i,]
        return (eval(fn, e, enclos))
    }

    if (dat$ns > 0) {
        firstval <- evalfor(1)
        res <- rep(firstval, dat$ns)
        if (dat$ns > 1) {
            for (i in 2:dat$ns) {
                res[i] <- evalfor(i)
            }
        }
        names(res) <- rownames(dat$x)
        return (res)
    } else {
        return (c())
    }
}

pcl.apply.f <- function(dat, expr, fn = substitute(expr), enclos = parent.frame()) {
    # Function application over features

    # Todo: add metaf support

    evalfor <- function(i) {
        e <- c(list(Name = colnames(dat$x)[i],
                    x = dat$x[,i],
                    Index = i),
               as.list(dat$meta),
               as.list(dat$x[,i]))
        return (eval(fn, e, enclos))
    }

    if (dat$nf > 0) {
        firstval <- evalfor(1)
        res <- rep(firstval, dat$nf)
        if (dat$nf > 1) {
            for (i in 2:dat$nf) {
                res[i] <- evalfor(i)
            }
        }
        names(res) <- colnames(dat$x)
        return (res)
    } else {
        return (c())
    }
}


# -----------------------------------------------------------------------------
#  PCL Filtering

pcl.only <- function(x, clade, rank, keepname=F) {
    # Slice out a subset of the features by clade and rank
    # Assumes features represent bug abundances stored in k__Bacteria format

    lx <- is.list(x)
    if (lx) {
        stopifnot(!("metaf" %in% names(x))) # This function has not yet had metaf support added
        l <- x
        x <- l$x
    }

    if (!missing(clade)) {
        # Pull out the columns with only that clade
        keep <- grepl(sprintf("(^|\\|)%s($|\\|)", clade), colnames(x))
        x <- x[,keep,drop=F]
        killpat <- sprintf("^(|.*\\|)%s", clade)
        if (!keepname) {
            colnames(x) <- mapply(function(n){gsub(killpat, clade, n)}, colnames(x))
        }
    }

    if (!missing(rank)) {
        # Pull out the columns with only that taxonomic rank
        keep <- grepl(sprintf("(^|\\|)%s__[^\\|]+$", rank), colnames(x))
        x <- x[,keep,drop=F]
        if (!keepname) {
            colnames(x) <- mapply(function(n){gsub("^.*\\|", "", n)}, colnames(x))
        }
    }

    if (lx) {
        l$x <- x
        l$nf <- dim(x)[2]
        return (l)
    } else {
        return (x)
    }
}

pcl.filter.s <- function(dat, keepwhat, keep, enclos = parent.frame()) {
    # Filter samples based on some criterion
    if (missing(keep)) {
        keep <- pcl.apply.s(dat, fn=substitute(keepwhat), enclos=enclos)
    }
    dat$x <- dat$x[keep,,drop=F]
    dat$meta <- dat$meta[keep,,drop=F]
    dat$ns <- dim(dat$x)[1]
    return (dat)
}

pcl.filter.f <- function(dat, keepwhat, keep, enclos = parent.frame()) {
    # Filter features based on some criterion
    if (missing(keep)) {
        keep <- pcl.apply.f(dat, fn=substitute(keepwhat), enclos=enclos)
    }
    dat$x <- dat$x[,keep,drop=F]
    if ("metaf" %in% names(dat)) {
        dat$metaf <- dat$metaf[keep,,drop=F]
    }
    dat$nf <- dim(dat$x)[2]
    return (dat)
}

pcl.top.s <- function(dat, expr=mean(x), n=10, fn=substitute(expr),
                      enclos = parent.frame(),
                      statistic=pcl.apply.s(dat, fn=fn, enclos=enclos)) {
    # Keep only the top n of the samples
    nonNA <- sum(!is.na(statistic))
    thresh <- if (nonNA == 0) { 0 } else {
        quantile(statistic, max(nonNA-n, 0)/nonNA, na.rm=T)
    }
    return (pcl.filter.s(dat, keep=!is.na(statistic) & statistic>=thresh))
}

pcl.top.f <- function(dat, expr=mean(x), n=10, fn=substitute(expr),
                      enclos = parent.frame(),
                      statistic=pcl.apply.f(dat, fn=fn, enclos=enclos)) {
    # Keep only the top n of the features
    nonNA <- sum(!is.na(statistic))
    thresh <- if (nonNA == 0) { 0 } else {
        quantile(statistic, max(nonNA-n, 0)/nonNA, na.rm=T)
    }
    return (pcl.filter.f(dat, keep=!is.na(statistic) & statistic>=thresh))
}

pcl.reorder.s <- function(dat, ord) {
    # Reorders the samples according to the new ordering
    ord <- ord[!is.na(ord)]
    dat$x <- dat$x[ord,,drop=F]
    dat$meta <- dat$meta[ord,,drop=F]
    dat$ns <- length(ord)
    return (dat)
}

pcl.reorder.f <- function(dat, ord, na.val=NA, keep.na=!is.na(na.val)) {
    # Reorders the features according to the new ordering
    setNames <- F
    if (!is.numeric(ord)) {
        featNames <- ord
        setNames <- T
        ord <- match(ord, colnames(dat$x))
    }
    if (!keep.na) {
        ord <- ord[!is.na(ord)]
    }
    dat$x <- dat$x[,ord,drop=F]
    if (keep.na && !is.na(na.val)) {
        dat$x[,is.na(ord)] <- na.val
    }
    if ("metaf" %in% names(dat)) {
        dat$metaf <- dat$metaf[ord,,drop=F]
    }
    dat$nf <- length(ord)
    if (setNames) {
        colnames(dat$x) <- featNames
    }
    return (dat)
}

pcl.match.f <- function(dat1, dat2, unmatched.val=NULL) {
    # Matches dat1 to the columns of dat2

    stopifnot(!("metaf" %in% names(dat1))) # This function has not yet had metaf support added

    mt <- match(colnames(dat2$x), colnames(dat1$x))
    if (!is.null(unmatched.val)) {
        datm <- pcl.reorder.f(dat1, mt, keep.na=T)
        datm$x[,is.na(mt)] <- unmatched.val
        colnames(datm$x) <- colnames(dat2$x)
        return (datm)
    } else {
        return (pcl.reorder.f(dat1, mt))
    }
}

pcl.match.s <- function(dat1, dat2) {
    # Matches dat1 to the rows of dat2
    return (pcl.reorder.s(dat1, match(rownames(dat2$x), rownames(dat1$x))))
}

pcl.match <- function(dat1, dat2, unmatched.val=NULL) {
    # Matches dat1 to the rows and columns of dat2
    return (pcl.match.f(pcl.match.s(dat1, dat2), dat2, unmatched.val))
}

pcl.merge.meta <- function(dat, datm) {
    # Merges metadata from datm into dat
    # Data tables are untouched

    stopifnot(!("metaf" %in% names(datm))) # This function has not yet had metaf support added

    sm <- match(rownames(dat$meta), rownames(datm$meta))
    unm2m <- !(colnames(datm$meta) %in% colnames(dat$meta))
    dat$meta <- cbind(dat$meta, datm$meta[sm, which(unm2m), drop=F])

    return (dat)
}

pcl.merge <- function(..., lst=NULL) {
    # Merges a bunch of pcl's

    merge2 <- function(dat1, dat2) {
        # Merges two pcl's

        # This function has not yet had metaf support added
        stopifnot(!(("metaf" %in% names(dat1)) || ("metaf" %in% names(dat2))))

        # Find the features, metadata and samples that are common
        fm <- match(colnames(dat1$x), colnames(dat2$x))
        mm <- match(colnames(dat1$meta), colnames(dat2$meta))
        sm <- match(rownames(dat1$x), rownames(dat2$x))

        # Find the features, metadata and samples in dat2 that aren't in dat1
        unm2f <- !(colnames(dat2$x) %in% colnames(dat1$x))
        unm2m <- !(colnames(dat2$meta) %in% colnames(dat1$meta))
        unm2s <- !(rownames(dat2$x) %in% rownames(dat1$x))

        if (any(!is.na(fm)) && any(!is.na(sm))) {
            # There are samples AND features in common
            # Make sure they're equal
            stopifnot(all(dat1$x[!is.na(sm), !is.na(fm), drop=F] ==
                          dat2$x[sm[!is.na(sm)], fm[!is.na(fm)], drop=F]),
                      "Overlapping data must be equal to merge PCL's")
        }

        if (any(!is.na(mm)) && any(!is.na(sm))) {
            # There are samples AND metadata in common
            # Make sure they're equal
            stopifnot(all(dat1$meta[!is.na(sm), !is.na(mm), drop=F] ==
                          dat2$meta[sm[!is.na(sm)], mm[!is.na(mm)], drop=F]))
        }

        # All overlapping values are consistent between the two tables - merge!

        # Create the new merged dat
        dat <- list(ns = dat1$ns + sum(unm2s),
                    nf = dat1$nf + sum(unm2f))

        # Merge the data tables
        dat$x <- rbind(
            cbind(dat1$x, dat2$x[sm, which(unm2f), drop=F]),
            dat2$x[which(unm2s), c(fm, which(unm2f)), drop=F])

        # Merge the metadata
        # NOTE: apparently data frames cannot be indexed with NAs, so the statement
        #       above for the data matrix cannot be used, and the bottom left part
        #       of the metadata data frame must be generated separately:
        bl <- dat1$meta[rep(1,sum(unm2s)),, drop=F]
        bl[,is.na(mm)] <- NA
        bl[,!is.na(mm)] <- dat2$meta[which(unm2s), mm[!is.na(mm)], drop=F]
        rownames(bl) <- rownames(dat2$meta)[which(unm2s)]
        if (dim(dat1$meta)[2] == 0) {
            # NOTE 2: rbind with a 0-column data.matrix does not produce the
            #         expected result.. hence this special case:
            dat$meta <- dat2$meta[c(sm, which(unm2s)), which(unm2m), drop=F]
        } else {
            dat$meta <- cbind(rbind(dat1$meta, bl),
                dat2$meta[c(sm, which(unm2s)), which(unm2m), drop=F])
        }

        return (dat)
    }

    lst <- c(list(...), lst)
    if(length(lst) == 0) {
        return (pcl.empty)
    } else {
        merge_sublist <- function(lst) {
            if (length(lst) > 1) {
                midpt <- floor(length(lst)/2)
                return (merge2(merge_sublist(lst[1:midpt]),
                               merge_sublist(lst[(midpt+1):length(lst)])))
            } else {
                stopifnot(length(lst) == 1)
                return (lst[[1]])
            }
        }
        return (merge_sublist(lst))
    }
}

pcl.group.f <- function(dat, map, avg=F, mapping=F) {
    # Group features according to a mapping
    # If avg is T, features are averaged rather than summed
    # If mapping is T, then map is treated as a named vector providing
    # a feature name -> output name map, rather than a direct output name for
    # each feature.

    library(dplyr)

    if (mapping) {
        map <- map[match(colnames(dat$x), names(map))]
    }
    xm <- split(as.data.frame(t(dat$x)), map) %>%
        lapply(colwise(if (avg) {mean} else {sum})) %>%
        do.call(rbind, .) %>% t

    if ("metaf" %in% names(dat)) {
        allsame <- function(x) length(x)==0 || all(x==x[1])
        mfm <- split(pcl$metaf, map) %>%
            lapply(colwise(function(x)
                if (any(!is.na(x)) && allsame(x[!is.na(x)]))
                    {x[min(which(!is.na(x)))]} else {NA})) %>%
            do.call(rbind, .)

        dat$metaf <- mfm
    }

    dat$x <- xm
    dat$nf <- ncol(xm)

    return (dat)
}

pcl.group.s <- function(dat, map, avg=F) {
    # Group samples according to a mapping
    # If avg is T, samples are averaged rather than summed
    # Metadata which is consistent across all grouped samples is kept
    # inconsistent metadata is tossed

    library(dplyr)

    xm <- split(as.data.frame(dat$x), map) %>%
        lapply(colwise(if (avg) {mean} else {sum})) %>%
        do.call(rbind, .)

    metaclean <- dat$meta[,colnames(dat$meta) != ""]
    mm <- split(metaclean, map) %>%
        lapply(., colwise(function(x)
            if (any(!is.na(x)) && all(x==na.omit(x)[1], na.rm=T))
            {na.omit(x)[1]} else {NA})) %>%
        do.call(rbind, .)

    dat$x <- xm
    dat$meta <- mm
    dat$ns <- nrow(xm)

    return (dat)
}

pcl.t <- function(dat) {
    return (pcl.make(t(dat$x)))
}

# -----------------------------------------------------------------------------
#  Misc PCL Helpers

pcl.map.fnames <- function(dat, mapexpr, enclos = parent.frame()) {
    # Map a function to the feature names
    colnames(dat$x) <- pcl.apply.f(dat, fn=substitute(mapexpr), enclos=enclos)
    return (dat)
}

pcl.map.snames <- function(dat, mapexprenclos = parent.frame()) {
    # Map a function to the sample names
    rownames(dat$x) <- pcl.apply.s(dat, fn=substitute(mapexpr), enclos=enclos)
    if (dim(dat$meta)[2] > 0) {
        rownames(dat$meta) <- rownames(dat$x)
    } else {
        dat$meta <- data.frame(row.names = rownames(dat$x))
    }
    return (dat)
}

pcl.nicenames <- function(dat) {
    # Prettifies feature names by removing stratification and transforming
    # s__Xyz_Abc_1 to Xyz Abc 1
    # Works on pathway names as well as taxonomic trees
    # Also works on ;-separated taxonomies
    return (pcl.map.fnames(dat, gsub("_", " ", gsub(
        "^((.*\\|)?\\w__(.*\\.\\w__)?)|^.*\\||^.*; ", "", Name))))
}

pcl.normalize <- function(dat, s=rowSums(dat$x, na.rm=T)) {
    # Normalize samples to sum to 1
    s[s==0] <- 1
    if (dat$nf >= 1)
        dat$x <- dat$x / s
    return (dat)
}

pcl.sort.f <- function(dat, feat, fn=substitute(feat), enclos = parent.frame()) {
    return (pcl.reorder.f(dat, order(pcl.apply.f(dat, fn=fn, enclos=enclos))))
}

pcl.sort.s <- function(dat, feat, fn=substitute(feat), enclos = parent.frame()) {
    return (pcl.reorder.s(dat, order(pcl.apply.s(dat, fn=fn, enclos=enclos))))
}


# -----------------------------------------------------------------------------
#  Common analysis functions

pcl.compmatrix.s <- function(dat, fun, symmetric=TRUE) {
    # Sample-sample comparison matrix
    # FUNCTION LARGELY UNTESTED..

    val1 <- fun(dat$x[1,], dat$x[1,]);
    comp <- matrix(val1, dat$ns, dat$ns)
    for (j in 2:dat$ns) {
        comp[1,j] <- fun(dat$x[1,], dat$x[j,])
    }
    if (symmetric) {
        for (i in 2:dat$ns) {
            for (j in i:dat$ns) {
                comp[i,j] <- fun(dat$x[i,], dat$x[j,])
                comp[j,i] <- comp[i,j]
            }
        }
    } else {
        for (i in 2:dat$ns) {
            for (j in 1:dat$ns) {
                comp[i,j] <- fun(dat$x[i,], dat$x[j,])
            }
        }
    }
    return (comp)
}

# -----------------------------------------------------------------------------
#  Heatmap visualization function

pcl.heatmap <- function(dat, ..., meta = F, logspace = F, pseudocount = 0,
                        zerocolor = NA, sqrtspace = F, gamma = 2, as = F,
                        ev = 10/(dat$ns*dat$nf), # Extreme values
                        minx = quantile(dat$x[is.finite(dat$x) & (dat$x>0)], ev),
                        maxx = quantile(dat$x[is.finite(dat$x) & (dat$x>0)], 1-ev),
                        metanames = NA, reorder = T, divergent0 = NA) {

    if (dat$ns == 0 || dat$nf == 0) {
        return (plot.new())
    }

    # Pretty heatmap of the data in dat
    library(pheatmap)

    params <- list(...)

    neg <- any(dat$x[is.finite(dat$x)]<0)
    if (is.na(divergent0)) {
        divergent0 <- neg
    }

    def.params <- list()
    def.params$show_colnames <- dat$ns <= 50
    def.params$show_rownames <- dat$nf <= 40
    def.params$clustering_distance_rows <- if(divergent0) {"euclidean"} else {"bray/curtis"}
    def.params$clustering_distance_cols <- def.params$clustering_distance_rows
    def.params$cluster_rows <- dat$nf > 1
    def.params$cluster_cols <- dat$ns > 1
    #def.params$color <- colorRampPalette(c("black","blue","cyan","yellow","red"))(100)
    if(divergent0) {
        def.params$color <- colorRampPalette(c("dodgerblue","white","firebrick"))(101)
        mv <- max(1e-10, max(-max(minx, min(dat$x, na.rm=T)), min(maxx, max(dat$x, na.rm=T))))
        def.params$breaks <- c(seq(from=-mv, to=mv, length=101))
    } else {
        library(viridis)
        def.params$color <- viridis(100)
        #def.params$color <- colorRampPalette(c("dodgerblue","goldenrod1","firebrick"))(100)
    }
    def.params$treeheight_row <- 100
    def.params$treeheight_col <- 100
    def.params$border_color <- NA
    def.params$cuttree_rows <- 5
    def.params$cuttree_cols <- 5

    x <- t(dat$x)
    x[!is.finite(x)] <- 0
    zero <- x == 0

    # Override the clustering
    override_clustering <- function(mat, distance, def) {
        if(class(distance) == "dist") {
            return (distance)
        } else {
            method <- if(is.null(distance)){def}else{distance}
            if (as) {
                mat <- asin(sqrt(mat/100))
            }
            if (method == "correlation") {
                D <- as.dist(1 - cor(t(mat)))
            } else if (method %in% c("manhattan", "euclidean", "canberra",
                                     "bray", "kulczynski", "jaccard", "gower",
                                     "altGower", "morisita", "horn", "mountford",
                                     "raup" , "binomial", "chao", "cao", "mahalanobis")) {
                library(vegan)
                D <- vegdist(mat, method)
            } else if (method %in% c("steinhaus", "sorensen", "ochiai",
                                     "ruzicka", "bray/curtis", "roberts", "chisq")) {
                library(labdsv)
                D <- dsvdis(mat, method)
            } else {
                D <- dist(mat, method = method)
            }
            return (D)
        }
    }

    params$clustering_distance_rows <- override_clustering(x,
        params$clustering_distance_rows, def.params$clustering_distance_rows)
    params$clustering_distance_cols <- override_clustering(t(x),
        params$clustering_distance_cols, def.params$clustering_distance_cols)

    # Perform data transformations
    x <- pmin(pmax(x, minx), maxx) + pseudocount
    if (logspace) {
        x <- log10(x)
    } else if (sqrtspace) {
        x <- x ** (1/gamma)
    }

    # Make zero special
    if (!divergent0 && (length(zerocolor)>1 || !is.na(zerocolor) && any(zero))) {
        x[zero] <- min(x) - (max(x)-min(x)) * (2/length(def.params$color))
        if ("color" %in% names(params)) {
            params$color <- c(if (is.character(zerocolor)){zerocolor}else{"gray95"},
                              params$color[1], params$color)
        } else {
            def.params$color <- c(if (is.character(zerocolor)){zerocolor}else{"gray95"},
                                  def.params$color[1], def.params$color)
        }
    }

    # Add metadata annotation
    if (length(meta) > 0 && (length(meta) > 1 || !is.na(meta))) {
        if (length(meta) == 1 && is.logical(meta)) {
            if (meta) {
                def.params$annotation_col <- dat$meta[,, drop = F]
            }
        } else {
            def.params$annotation_col <- dat$meta[,match(meta, colnames(dat$meta)), drop = F]
        }

        for (i in seq_along(def.params$annotation_col)) {
            if (class(def.params$annotation_col[,i]) == "logical") {
                def.params$annotation_col[,i] <- factor(def.params$annotation_col[,i])
            }
        }

        if (length(metanames) > 1 || !is.na(metanames)) {
            names2 <- names(def.params$annotation_col)
            nmatch <- match(names2, names(metanames))
            names2[!is.na(nmatch)] <- metanames[nmatch[!is.na(nmatch)]]
            names(def.params$annotation_col) <- names2
        }
    }

    if (reorder) {
        def.params$clustering_callback <- function(hc, u) {
            library(seriation)
            m <- if (all(dim(u) == dim(x)) && all(u == x)) {
                params$clustering_distance_rows
            } else {
                params$clustering_distance_cols
            }
            ord <- get_order(seriate(m, control=list(hclust=hc), method="OLO"), 1)
            hc$order <- ord
            return (hc)
        }
    }

#     def.params$clustering_callback <- function(hc, x) {
#         library(vegan)
#         D <- dsvdis(x, "bray/curtis")
#         return (reorder.hclust(hc, D)) # D is not what is expected here
#     }

    args <- c(list(mat=x), params, def.params)
    args <- args[unique(names(args))]
    do.call(pheatmap, args)
}

# -----------------------------------------------------------------------------
#  Ordination functions

pcl.nmds <- function(dat, asinsqrt = T, remove_zeros = T, ...) {
    # Performs a NMDS ordination

    x <- dat$x

    # Arcsin Sqrt transform
    if (asinsqrt) {
        maxx <- max(x)
        if (maxx > 1) {
            if (maxx > 100) {
                x <- x / maxx
            } else {
                x <- x / 100
            }
        }
        x <- asin(sqrt(x))
    }

    # Warn on zero samples to avoid some confusion about the wierd error
    # vegan gives in this condition
    zerosamples <- apply(x, 1, function(xx){all(xx==0)})
    if (any(zerosamples)) {
        if (remove_zeros) {
            x <- x[!zerosamples,]
        } else {
            warning("There are samples with zero abundances. This may cause errors in metaMDS.")
        }
    }

    library(vegan)
    ord <- metaMDS(x, autotransform = F, ...)

    return (ord)
}

pcl.pcoa <- function(dat, D = NA, as = F, asinsqrt = as, index = "bray/curtis", k = 2) {
    # Performs a PCoA ordination

    library(labdsv)

    if (length(D) == 1 && is.na(D)) {
        x <- dat$x
        x[!is.finite(x)] <- 0

        # Arcsin Sqrt transform
        if (asinsqrt) {
            maxx <- max(x)
            if (maxx > 1) {
                if (maxx > 100) {
                    x <- x / maxx
                } else {
                    x <- x / 100
                }
            }
            x <- asin(sqrt(x))
        }

        k <- max(k, 2)
        D <- dsvdis(x, index)
    }

    pc <- pco(D, k)

    ord <- pc
    toteig <- sum(pc$eig[pc$eig>0])
    ord$ordnames <- sprintf("PCo %d (%0.1f%%)", 1:k, 100 * pc$eig[1:k] / toteig)
    return (ord)
}

pcl.pcoa_sharpel <- function(dat, D = NA, as = F, asinsqrt = as, index = "bray/curtis", k = 2, var_capture=0.9) {
    # Performs a PCoA ordination

    library(labdsv)

    if (length(D) == 1 && is.na(D)) {
        x <- dat$x
        x[!is.finite(x)] <- 0

        # Arcsin Sqrt transform
        if (asinsqrt) {
            maxx <- max(x)
            if (maxx > 1) {
                if (maxx > 100) {
                    x <- x / maxx
                } else {
                    x <- x / 100
                }
            }
            x <- asin(sqrt(x))
        }

        k <- max(k, 2)
        D <- dsvdis(x, index)
    }

    pc <- pco(D, dat$ns-1)

    varfrac <- cumsum(pmax(0, pc$eig))
    varfrac <- varfrac / max(varfrac)
    qk <- min(which(varfrac > var_capture))

    library(pcaL1)
    spc <- sharpel1pca(sweep(pc$points[,1:qk], 2, varfrac[1:qk], FUN="*"), k, projections="l1")
    ord <- list(points=spc$scores,
                ordnames=sprintf("Axis %d (%0.1f%%)", 1:k, 100 * spc$dispExp[1:k] * varfrac[qk]))

    return (ord)
}

pcl.pca_sharpel <- function(dat, k = 2) {
    library(pcaL1)
    spc <- sharpel1pca(dat$x, k, projections="l1")
    ord <- list(points=spc$scores,
                ordnames=sprintf("Axis %d (%0.1f%%)", 1:k, 100 * spc$dispExp[1:k]))

    return (ord)
}

pcl.dmap <- function(dat, D = NA, as = F, asinsqrt = as, index = "bray/curtis", k = 2) {
    # Performs a PCoA ordination

    library(diffusionMap)
    library(labdsv)

    if (length(D) == 1 && is.na(D)) {
        x <- dat$x
        x[!is.finite(x)] <- 0

        # Arcsin Sqrt transform
        if (asinsqrt) {
            maxx <- max(x)
            if (maxx > 1) {
                if (maxx > 100) {
                    x <- x / maxx
                } else {
                    x <- x / 100
                }
            }
            x <- asin(sqrt(x))
        }

        k <- max(k, 2)
        D <- dsvdis(x, index)
    }

    dm <- diffuse(D, maxdim=k)

    ord <- dm
    ord$points <- dm$X
    ord$ordnames <- sprintf("DMap %d", 1:k)
    rownames(ord$points) <- rownames(dat$x)
    return (ord)
}

pcl.lda <- function(dat, meta, as = T, asinsqrt = as, k = 2) {
    # Linear Discriminant Analysis

    x <- dat$x
    groups <- dat$meta[,meta]

    if (any(is.na(groups))) {
        warning("NAs in the metadata")
        x <- x[!is.na(groups),]
        groups <- groups[!is.na(groups)]
    }

    # Arcsin Sqrt transform
    if (asinsqrt) {
        maxx <- max(x)
        if (maxx > 1) {
            if (maxx > 100) {
                x <- x / maxx
            } else {
                x <- x / 100
            }
        }
        x <- asin(sqrt(x))
    }

    # lda can't handle variables with zero variance.. remove them
    for (i in seq_along(levels(groups))) {
        xvar <- apply(x[groups==levels(groups)[i],], 2, var)
        x <- x[,xvar > 0]
    }

    library(MASS)
    res <- lda(x, groups)
    p <- predict(res, x)

    ord <- list(lda_res = res)
    if (dim(p$x)[2] >= 2) {
        ord$points <- p$x
        ord$ordnames <- sprintf("LD %d", 1:dim(p$x)[2])
    } else {
        ord2 <- pcl.pcoa(dat, k=1)

        ord$points <- cbind(p$x, ord2$points[,1])
        ord$ordnames <- c("LD 1", ord2$ordnames[1])
    }

    return (ord)
}

pcl.tsne <- function(dat, D = NA, as = F, asinsqrt = as,index = "bray/curtis", k = 2, seed = 1234) {
    # Performs a PCoA ordination

    library(tsne)
    library(labdsv)

    if (length(D) == 1 && is.na(D)) {
        x <- dat$x
        x[!is.finite(x)] <- 0

        # Arcsin Sqrt transform
        if (asinsqrt) {
            maxx <- max(x)
            if (maxx > 1) {
                if (maxx > 100) {
                    x <- x / maxx
                } else {
                    x <- x / 100
                }
            }
            x <- asin(sqrt(x))
        }

        k <- max(k, 2)
        D <- dsvdis(x, index)
    }

    k <- max(k, 2)

    set.seed(seed)
    tsn <- tsne(D, k=k)
    ord <- list(points = tsn, ordnames = sprintf("t-SNE %d", 1:k))
    rownames(ord$points) <- rownames(dat$x)
    return (ord)
}

pcl.ordplot <- function(dat, ord, pcos = 2, pointoutline = T,
        colour = NA, colour_title = NA, colour_names = NA, colour_override = NA,
        shape = NA,  shape_title = NA,  shape_names = NA,  shape_override = NA,
        size = NA, size_title = NA, size_names = NA, size_override = NA,
        size_abs = NA, outline_size = NA,
        surface = NA, surf_maj_lines = 8, surf_min_lines = 4, surf_extent = 1,
        surf_smoothness = 1, surf_quality = 100,
        centroid = NA, loading = NA, sigloadings = NA,
        enriched = NA, text_halo = T,
        colour_log = F, colour_log_pseudocount = 0, colour_log_base = 10,
        arrows_fixed = NULL, arrows = NULL, arrow_text = T, arrow_sqrtnorm = T,
        connect = NA, sequence = NA, connectwidth = NA, sortby = NA, decreasing = F) {

    # ggplot2-based function to visualize the results of an ordination
    library(ggplot2)

    # Set up the plot
    ggp <- ggplot() + theme_classic() +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    gptopt <- list()
    gptaes <- list(x="dim1", y="dim2")

    # Get the dimensions to plot, and set the axes appropriately
    if (length(pcos) == 1) {
        pcos <- max(pcos,2)
        pcos <- c(pcos-1, pcos)
    } else {
        stopifnot(length(pcos) == 2)
    }
    stopifnot(max(pcos) <= dim(ord$points)[2])
    pts <- as.data.frame(ord$points[,pcos])
    if (is.null(ord$ordnames)) {
        dims <- colnames(pts)
    } else {
        dims <- ord$ordnames[pcos]
    }
    colnames(pts) <- c("dim1", "dim2")
    ggp <- ggp + xlab(dims[1]) + ylab(dims[2])

    # Use this mapping to allow the ordination to be done on a subset of the
    # data without requiring that dat also be that subset
    metamatch <- match(rownames(pts), rownames(dat$x))
    if (any(is.na(metamatch))) {
        pts <- pts[!is.na(metamatch),]
        metamatch <- metamatch[!is.na(metamatch)]
    }

    getmeta <- function(name) {
        tgt <- name == colnames(dat$meta)
        if (any(tgt)) {
            return (dat$meta[metamatch,tgt])
        }
        tgt <- name == colnames(dat$x)
        if (any(tgt)) {
            return (dat$x[metamatch,min(which(tgt))])
        }

        stop(sprintf("%s is not a metadata or feature.", name))
    }

    # Apply sorting
    if (!is.na(sortby)) {
        s <- getmeta(sortby)
        pts <- pts[order(s, decreasing=decreasing),]
        metamatch <- match(rownames(pts), rownames(dat$x))
    }

    # Draw connectors first
    if (!is.na(connect) || !is.na(sequence)) {
        df <- data.frame(x=pts$dim1, y=pts$dim2)

        if (!is.na(connect)) {
            connectby <- getmeta(connect)
            df$group <- connectby
        }
        if (!is.na(sequence)) {
            seqby <- getmeta(sequence)
            df <- df[order(seqby),]
        }
        if (!is.na(connect)) {
            ggp <- ggp + geom_path(data=df, aes(x=x, y=y, group=group), size=ifelse(is.na(connectwidth), 1, connectwidth))
        } else {
            ggp <- ggp + geom_path(data=df, aes(x=x, y=y), size=ifelse(is.na(connectwidth), 1, connectwidth))
        }
    }

    #     set_by <- function(bywhat, names, override, title,
    #                        aeskey, default, ggscale_manual, legoverride) {
    #         if (!is.na(bywhat)) {
    #             byval <- dat$meta[[bywhat]][metamatch]
    #             if (length(names) > 1 || !is.na(names)) {
    #                 # reorder the factors
    #                 byval <- factor(levels(shapeby)[as.numeric(shapeby)],
    #                                   levels=names(names))
    #             }
    #             pts[[aeskey]] <- byval
    #             gptaes[[aeskey]] <- aeskey
    #             if (is.na(title)) {
    #                 title <- bywhat
    #             }
    #
    #             guideopt <- list()
    #             guideopt[[aeskey]] <- guide_legend(title = shape_title,
    #                                                override.aes = legoverride)
    #             ggp <- ggp + do.call(guides, guideopt)
    #             if (length(override) > 1 || !is.na(override)) {
    #                 ggp <- ggp + ggscale_manual(values=unlist(override),
    #                                             breaks=names(names),
    #                                             labels=unlist(names))
    #             }
    #
    #         } else if (length(override) == 1 && !is.na(override)) {
    #             gptopt[[aeskey]] <- override
    #         } else {
    #             gptopt[[aeskey]] <- default
    #         }
    #     }

    # Draw the surface
    if (!is.na(surface)) {
        surfby <- getmeta(surface)

        # pick an appropriate bandwidth
        bw1 <- surf_smoothness * 1.06 * sd(pts$dim1) * length(pts$dim1)^-0.2
        bw2 <- surf_smoothness * 1.06 * sd(pts$dim2) * length(pts$dim2)^-0.2

        library(signal) # for unwrap

        # only make the surface from finite values
        keep <- !is.na(surfby) & is.finite(surfby)
        surfby <- surfby[keep]
        xp <- pts$dim1[keep]
        yp <- pts$dim2[keep]

        # distance from the points to draw
        extent <- 0.5 * surf_extent^2 / surf_smoothness

        # surface from a Gaussian weighted mean
        gkmeansurf <- function(x, y) {
            dx <- x - xp
            dy <- y - yp
            chi2 <- (dx / bw1)^2 + (dy / bw2)^2
            kern <- exp(-chi2)

            theta <- unwrap(atan2(dy, dx))
            mnchi2 <- min(chi2)
            if (mnchi2 > extent && (max(theta) - min(theta)) < pi) { return (NA) }

            return (-sum(surfby * kern) / sum(kern))
        }

        # drop a grid on the ordination
        n <- surf_quality + 1
        border <- surf_extent
        xx <- matrix(seq(from=min(pts$dim1) - border * bw1,
                         to  =max(pts$dim1) + border * bw1,
                         length=n), n, n, byrow=F)
        yy <- matrix(seq(from=min(pts$dim2) - border * bw2,
                         to  =max(pts$dim2) + border * bw2,
                         length=n), n, n, byrow=T)

        # measure the surface at all points
        z <- mapply(gkmeansurf, as.vector(xx), as.vector(yy))

        # only keep points in range
        keep <- !is.na(z)

        # plot contours
        zrange <- max(z[keep]) - min(z[keep])
        cdata <- data.frame(x=as.vector(xx)[keep], y=as.vector(yy)[keep], z=z[keep])
        ggp <- ggp + stat_contour(aes(x=x, y=y, z=z), data=cdata,
               size=0.5, colour="grey50", binwidth=zrange/(surf_maj_lines*(surf_min_lines+1)))
        ggp <- ggp + stat_contour(aes(x=x, y=y, z=z), data=cdata,
               size=1, colour="black", binwidth=zrange/surf_maj_lines)
    }

    # Set up the plot to shape by some metadata
    if (!is.na(shape)) {
        shapeby <- getmeta(shape)
        if (length(shape_names) > 1 || !is.na(shape_names)) {
            # reorder the factors to match the order in shape_names
            # so that the legend is displayed in the correct order
            shapeby <- factor(levels(shapeby)[as.numeric(shapeby)],
                              levels=names(shape_names))
        }
        pts$shape <- shapeby
        gptaes$shape <- "shape"
        if (is.na(shape_title)) {
            shape_title <- shape
        }

        ggp <- ggp + guides(shape = guide_legend(title = shape_title,
                                                 override.aes = list(size=4, fill="gray")))

        if (length(shape_override) > 1 || !is.na(shape_override)) {
            # User-defined shapes
            if (length(shape_names) == 1 && is.na(shape_names)) {
                shape_names <- shape_override
                shape_names[names(shape_override)] <- names(shape_override)
            }
            ggp <- ggp + scale_shape_manual(values=unlist(shape_override),
                                            breaks=names(shape_names),
                                            labels=unlist(shape_names))
        } else if (pointoutline) {
            # Automatic shapes - select the "fillable" ones
            shps <- list()
            shapeby <- factor(shapeby)
            for (i in 1:length(levels(shapeby))) {
                shps[[levels(shapeby)[i]]] <- 20+i
            }
            ggp <- ggp + scale_shape_manual(values=unlist(shps))
        } else {
            # Automatic shapes - no outlines
            shps <- list()
            shapeby <- factor(shapeby)
            for (i in 1:length(levels(shapeby))) {
                shps[[levels(shapeby)[i]]] <- 15+i
            }
            ggp <- ggp + scale_shape_manual(values=unlist(shps))
        }

    } else if (length(shape_override) == 1 && !is.na(shape_override)) {
        gptopt$shape <- shape_override
    } else {
        gptopt$shape <- if (pointoutline) {21} else {19}
    }

    # Set up the plot to size by some metadata
    if (!is.na(size)) {
        sizeby <- getmeta(size)

        gptaes$size <- "size"
        if (is.na(size_title)) {
            size_title <- size
        }

        if (is.numeric(sizeby)) {
            if (!is.na(size_abs)) {
                sizeby <- sizeby * size_abs
            }
        } else if (length(size_override) > 1 || !is.na(size_override)) {
            # User-defined shapes
            if (length(size_names) == 1 && is.na(size_names)) {
                size_names <- size_override
                size_names[names(size_override)] <- names(size_override)
            }
            ggp <- ggp + scale_size_manual(values=unlist(size_override),
                                            breaks=names(size_names),
                                            labels=unlist(size_names))
        }

        pts$size <- sizeby
        ggp <- ggp + guides(size = guide_legend(title = size_title))
    } else {
        if (is.na(size_abs)) {
            gptopt$size <- min(6, 50 / sqrt(length(metamatch)))
        } else {
            gptopt$size <- size_abs
        }
    }

    # Set up the plot to colour based on metadata
    if (pointoutline) {
        ggscale_manual <- scale_fill_manual
        ggscale_gradientn <- scale_fill_gradientn
        ggscale_brewer <- scale_fill_brewer
        ggcolorfield <- "fill"
        if (!is.na(outline_size)) {
            gptopt[["stroke"]] <- outline_size
        }
    } else {
        ggscale_manual <- scale_colour_manual
        ggscale_gradientn <- scale_colour_gradientn
        ggscale_brewer <- scale_colour_brewer
        ggcolorfield <- "colour"
    }
    if (!is.na(colour)) {
        colby <- getmeta(colour)
        if (length(colour_names) > 1 || !is.na(colour_names)) {
            # reorder the factors to match the order in colour_names
            # so that the legend is displayed in the correct order
            colby <- factor(as.character(colby),
                            levels=names(colour_names))
        }
        if (is.numeric(colby) && colour_log) {
            colby <- log(colby + colour_log_pseudocount, base=colour_log_base)
        }
        pts$colour <- colby
        gptaes[[ggcolorfield]] <- "colour"

        if (is.na(colour_title)) {
            colour_title <- colour
        }

        if (length(colour_override) > 1 || !is.na(colour_override)) {
            # Manual colouring
            if (length(colour_names) == 1 && is.na(colour_names)) {
                colour_names <- colour_override
                colour_names[names(colour_override)] <- names(colour_override)
            }
            if (is.numeric(colby)) {
                ggp <- ggp + ggscale_gradientn(colours = colour_override, na.value="grey80")
            } else {
                ggp <- ggp + ggscale_manual(values=unlist(colour_override),
                                            breaks=names(colour_names),
                                            labels=unlist(colour_names))
            }
            col_guide <- guide_legend(title = colour_title,
                                      override.aes = list(size=3, shape=21))
        } else if (is.numeric(colby)) {
            # Continuous data uses the RdYlGn palette
            #library(RColorBrewer)
            #ggp <- ggp + ggscale_gradientn(colours = colorRampPalette(c("dodgerblue","goldenrod1","firebrick"))(100))
            library(viridis)
            ggp <- ggp + ggscale_gradientn(colours = viridis(100), na.value="grey80")
                   #brewer.pal(n = 11, name = "YlGnBu"))
            #ggp <- ggp + ggscale_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlGn")))
            col_guide <- guide_colourbar(title = colour_title)
        } else if (length(levels(colby)) > 9) {
            library(RColorBrewer)
            ggp <- ggp + ggscale_manual(values = colorRampPalette(
                rev(brewer.pal(n = 7, name = "RdYlBu")))(length(levels(colby))))
            col_guide <- guide_legend(title = colour_title,
                                    override.aes = list(size=3, shape=21))
        } else {
            # Discrete data uses the Set1 palette
            ggp <- ggp + scale_fill_brewer(palette = "Set1")
            col_guide <- guide_legend(title = colour_title,
                                override.aes = list(size=3, shape=21))
        }
        if (pointoutline) {
            ggp <- ggp + guides(fill = col_guide)
        } else {
            ggp <- ggp + guides(colour = col_guide)
        }
    } else if (length(colour_override) == 1 && !is.na(colour_override)) {
        gptopt[[ggcolorfield]] <- colour_override
    } else {
        gptopt[[ggcolorfield]] <- "darkslategray4"
    }

    # Draw the points
    gptopt$data <- pts
    gptopt <- c(list(do.call(aes_string, gptaes)), gptopt)
    ggp <- ggp + do.call(geom_point, gptopt)

    # Helper to add text annotations
    xr <- 0.0018
    dx <- xr * (max(pts[,1]) - min(pts[,1]))
    dy <- 1.1 * xr * (max(pts[,2]) - min(pts[,2]))
    halo_quality <- 12
    put_text <- function(x, y, text, ...) {
        if (text_halo) {
            phi <- 2 * pi / halo_quality
            for (i in 1:halo_quality) {
                ggp <<- ggp + annotate("text", label=text, color="white",
                                       x = x + dx*cos((i - 0.5) * phi),
                                       y = y + dy*sin((i - 0.5) * phi), ...)
            }
        }
        ggp <<- ggp + annotate("text", label=text, x=x, y=y, ..., color="black")
    }

    # Add metadata centroids
    if (length(centroid) > 1 || !is.na(centroid)) {
        for (i in seq_along(centroid)) {
            meta <- getmeta(centroid[i])

            if (is.numeric(meta)) {
                d1 <- sum(pts$dim1 * meta) / sum(meta)
                d2 <- sum(pts$dim2 * meta) / sum(meta)
                put_text(d1, d2, centroid[i])
            } else {
                meta <- factor(meta) # remove unused
                for (j in seq_along(levels(meta))) {
                    d1 <- mean(pts$dim1[meta == levels(meta)[j]])
                    d2 <- mean(pts$dim2[meta == levels(meta)[j]])
                    put_text(d1, d2, levels(meta)[j])
                }
            }
        }
    }

    # Add metadata enrichment modes
    if (length(enriched) > 1 || !is.na(enriched)) {
        library(neldermead)
        for (i in seq_along(enriched)) {
            meta <- getmeta(enriched[i])

            if (!is.numeric(meta)) {
                stop(sprintf("%s must be numeric to find its enrichment mode"))
            }

            # select appropriate bandwidths
            bw1 <- 1.06 * sd(pts$dim1) * length(pts$dim1)^-0.2
            bw2 <- 1.06 * sd(pts$dim2) * length(pts$dim2)^-0.2

            # only use finite values
            keep <- !is.na(meta) & is.finite(meta)
            meta <- meta[keep]
            xp <- pts$dim1[keep]
            yp <- pts$dim2[keep]

            enr <- function(xy) {
                chi2 <- ((xy[1] - xp) / bw1)^2 + ((xy[2] - yp) / bw2)^2
                kern <- exp(-chi2)
                nsamps <- sum(kern)
                modulator <- nsamps / (nsamps + 1)
                return (modulator * -sum(meta * kern) / sum(kern))
            }

            bestxy <- 0
            bestnenr <- Inf
            for (x0 in seq(from=min(pts$dim1), to=max(pts$dim1), length=5)) {
                for (y0 in seq(from=min(pts$dim2), to=max(pts$dim2), length=5)) {
                    res <- fminsearch(enr, c(x0, y0))
                    if (neldermead.get(res, "fopt") < bestnenr) {
                        bestxy <- neldermead.get(res, "xopt")
                        bestnenr <- neldermead.get(res, "fopt")
                    }
                }
            }

            put_text(bestxy[1], bestxy[2], enriched[i])
        }
    }

    # Add arrows
    if (!is.null(arrows) || !is.null(arrows_fixed)) {
        if (is.null(arrows_fixed)) {
            arr_pcos <- matrix(0, 0, 2, dimnames=list(c(), dims))
        } else {
            arr_pcos <- arrows_fixed[,pcos]
        }
        if (!is.null(arrows)) {
            if (is.numeric(arrows)) {
                arrSet <- c()
                feats <- c(colnames(dat$x))
                for (i in seq_along(feats)) {
                    meta <- getmeta(feats[i])

                    if (is.numeric(meta)) {
                        d1 <- sum(pts$dim1 * meta) / sum(meta)
                        d2 <- sum(pts$dim2 * meta) / sum(meta)
                        arr <- matrix(c(d1, d2), 1, 2)
                        rownames(arr) <- feats[i]
                        arrSet <- rbind(arrSet, arr)
                    } else {
                        meta <- factor(meta) # remove unused
                        for (j in seq_along(levels(meta))) {
                            d1 <- mean(pts$dim1[meta == levels(meta)[j]])
                            d2 <- mean(pts$dim2[meta == levels(meta)[j]])
                            arr <- matrix(c(d1, d2), 1, 2)
                            rownames(arr) <- levels(meta)[j]
                            arrSet <- rbind(arrSet, arr)
                        }
                    }
                }
                arrSz <- arrSet[,1]^2 + arrSet[,2]^2
                I <- order(-arrSz)
                arr_pcos <- rbind(arr_pcos, arrSet[I[1:arrows],])
            } else {
                for (i in seq_along(arrows)) {
                    meta <- getmeta(arrows[i])

                    if (is.numeric(meta)) {
                        d1 <- sum(pts$dim1 * meta) / sum(meta)
                        d2 <- sum(pts$dim2 * meta) / sum(meta)
                        arr <- matrix(c(d1, d2), 1, 2)
                        rownames(arr) <- arrows[i]
                        arr_pcos <- rbind(arr_pcos, arr)
                    } else {
                        meta <- factor(meta) # remove unused
                        for (j in seq_along(levels(meta))) {
                            d1 <- mean(pts$dim1[meta == levels(meta)[j]])
                            d2 <- mean(pts$dim2[meta == levels(meta)[j]])
                            arr <- matrix(c(d1, d2), 1, 2)
                            rownames(arr) <- levels(meta)[j]
                            arr_pcos <- rbind(arr_pcos, arr)
                        }
                    }
                }
            }
        }

        ptsRS <- rowSums(pts[,1:2]^2)
        maxRadius <- sqrt(max(ptsRS[is.finite(ptsRS)]))
        maxArrRadius <- sqrt(max(rowSums(arr_pcos^2)))
        arr_pcos <- arr_pcos * (0.85 * maxRadius / maxArrRadius)
        if (arrow_sqrtnorm) {
            arr_len <- sqrt(rowSums(arr_pcos^2))
            maxarr_len <- max(arr_len)
            narr_len <- maxarr_len * sqrt(arr_len / maxarr_len)
            arr_pcos <- arr_pcos * (narr_len / arr_len)
        }
        colnames(arr_pcos) <- paste("dim", 1:length(pcos), sep="")

        ggp <- ggp + geom_segment(data=as.data.frame(arr_pcos),
                                  aes(x=0, y=0, xend=dim1, yend=dim2),
                                  colour="red", arrow=arrow())

        if (arrow_text) {
            for (i in seq_len(nrow(arr_pcos))) {
                rr <- arr_pcos[i,]
                hjust <- 0.5
                vjust <- 0.5
                if (abs(rr[1]) > abs(rr[2])) {
                    hjust <- if (rr[1] > 0) {0} else {1}
                } else {
                    vjust <- if (rr[2] > 0) {0} else {1}
                }
                tc <- rr * (1 + 0.015 * maxRadius / sqrt(sum(rr^2)))
                put_text(tc[1], tc[2], rownames(arr_pcos)[i],
                         hjust=hjust, vjust=vjust, color='red')
            }
        }
    }

    return (ggp)
}

# -----------------------------------------------------------------------------
#  Efficient massive list/vector creation

blmerge <- function(bl, v) {
    # Merge a list/vector v into a list-of-v's bl
    # The result is a new bl
    # To create a new, empty bl, use list()
    # To finalize the large list/vector, unlist the bl. Elements are guaranteed
    # to be in the same order as if they had been c'd together by this function

    i <- 1
    while (i <= length(bl) && !is.null(bl[[i]])) {
        v <- c(bl[[i]], v)
        bl[[i]] <- NULL
        i <- i + 1
    }

    bl[[i]] <- v
    return (bl)
}

# -----------------------------------------------------------------------------
#  Match function using binary search

bmatch <- function(needle, haystack, le=F, nomatch=NA) {
    # Searches for needles in the vector/list haystack, which must be sorted
    # Returns the indices or nomatch if it does not exist
    # If le is TRUE, then the indices of the highest values that are <=needle
    # are returned (which may be zero all of haystack is greater than needle)
    #
    # If there are ties in haystack, then the index of the last element of tie
    # is returned

    lo <- rep(1, length(needle))
    hi <- rep(length(haystack), length(needle))

    # maintain an active set of needles which we have not yet found
    active <- seq_along(lo)
    while (length(active) > 0) {
        # get the midpoints of each active needle
        mid <- floor((lo[active] + hi[active]) / 2)
        midv <- haystack[mid]

        # limit the range of each active needle based on whether it's higher
        # than the midpoint of its current range in the haystack
        gt <- midv > needle[active]
        lo[active[!gt]] <- mid[!gt] + 1
        hi[active[ gt]] <- mid[ gt] - 1

        # deactivate found needles
        active <- active[lo[active] <= hi[active]]
    }

    # hi contains the indices we are interested in - replace inexact matches
    # with nomatch
    if (!le) {
        found <- hi > 0
        hi[!found] <- nomatch
        hi[which(found)[haystack[hi[found]] != needle[found]]] <- nomatch
    }

    return (hi)
}

bmatch_sort <- function(x) {
    # Sorts x by its names so that it can be used in bmatch

    x[order(names(x))]
}

# -----------------------------------------------------------------------------
#  Plotting helpers

geom_boxplot_n <- function(offset = 0, ...) {
    # Put the number of observations over a boxplot

    # Adapted from:
    # http://stackoverflow.com/questions/28846348/add-number-of-observations-per-group-in-ggplot2-boxplot
    return (stat_summary(fun.data = function(x){
        return(c(y = max(x) + offset, label = length(x)))
    }, geom = "text", position = position_dodge(width = 0.75), vjust="bottom", ...))
}



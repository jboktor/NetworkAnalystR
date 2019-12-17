##################################################
## R script for NetworkAnalyst
## Description: functions only for single gene expression data
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# read tab delimited file
# stored in dataSet list object
# can have many classes, stored in meta.info
ReadTabExpressData <- function(fileName) {
  dataSet <- .readTabData(fileName)

  # rename data to data.orig
  int.mat <- dataSet$data
  dataSet$data <- NULL

  msg <- paste("a total of ", ncol(int.mat), " samples and ", nrow(int.mat), " features were found")

  # remove NA, null
  row.nas <- apply(is.na(int.mat) | is.null(int.mat), 1, sum)
  good.inx <- row.nas / ncol(int.mat) < 0.5
  if (sum(!good.inx) > 0) {
    int.mat <- int.mat[good.inx, ]
    msg <- c(msg, paste("removed ", sum(!good.inx), " features with over 50% missing values"))
  }
  # remove constant values
  filter.val <- apply(int.mat, 1, IQR, na.rm = T)
  good.inx2 <- filter.val > 0
  if (sum(!good.inx2) > 0) {
    int.mat <- int.mat[good.inx2, ]
    msg <- c(msg, paste("removed ", sum(!good.inx2), " features with constant values"))
  }

  if (nrow(int.mat) > 5000) {
    filter.val <- filter.val[good.inx2]
    rk <- rank(-filter.val, ties.method = "random")

    var.num <- nrow(int.mat)
    kept.num <- 0.95 * var.num
    int.mat <- int.mat[rk < kept.num, ]
    # msg <- c(msg, paste("removed 5% features with near-constant values"));
  }

  minVal <- min(int.mat, na.rm = T)
  na.inx <- is.na(int.mat)
  if (sum(na.inx) > 0) {
    int.mat[na.inx] <- minVal / 2
    # msg <- c(msg, "the remaining", sum(na.inx), "missing variables were replaced with data min");
  }
  current.msg <<- paste(msg, collapse = "; ")
  data.proc <- RemoveDuplicates(int.mat, "mean", quiet = T)

  # save processed data for download user option
  write.csv(data.proc, file = "data_processed.csv")
  saveRDS(data.proc, "data.proc.rds")

  dataSet <<- dataSet
  return(1)
}

GetSampleNumber <- function() {
  data.proc <- readRDS("data.proc.rds")
  return(ncol(data.proc))
}

GetMetaInfo <- function() {
  return(colnames(dataSet$meta.info))
}

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# idType: INVEX supported ID types
# lvlOpt: "NA" to keep original, other values will merge original ID to entrez gene IDs

# return the total matched gene number
# note: unmapped IDs will be retained as
# original label (i.e. intergenic regions) in further analysis

PerformDataAnnot <- function(org, dataType, idType, lvlOpt) {
  data.org <<- org
  SetInitLib(org)

  dataSet$type <- dataType
  dataSet$id.orig <- dataSet$id.current <- idType
  dataSet$annotated <- F
  # should not contain duplicates, however sanity check
  data.proc <- readRDS("data.proc.rds")
  dataSet$data.anot <- data.proc

  if (org != "NA" & idType != "NA") {
    feature.vec <- rownames(data.proc)
    anot.id <- doAnnotation(feature.vec, idType)

    # dataSet$annotation <- anot.id;
    saveRDS(anot.id, "annotation.rds")

    hit.inx <- !is.na(anot.id)
    matched.len <- sum(hit.inx)
    perct <- round(matched.len / length(feature.vec), 3) * 100
    thresh <- 0.1 # previous value of 0.25 is causing challenges
    # for datasets like Ppromelas with low annotation quality
    if (matched.len < length(feature.vec) * thresh) {
      current.msg <<- paste("Only ", perct, "% ID were matched. You may want to choose another ID type or use default.", sep = "")
    } else {
      current.msg <<- paste("ID annotation: ", "Total [", length(anot.id),
        "] Matched [", matched.len, "] Unmatched [", sum(!hit.inx), "]",
        collapse = "\n"
      )

      if (lvlOpt != "NA" | idType == "entrez") {
        # do actual summarization to gene level

        matched.entrez <- anot.id[hit.inx]
        data.anot <- data.proc[hit.inx, ]
        rownames(data.anot) <- matched.entrez
        current.msg <<- paste(current.msg, "Data is now transformed to gene-level (Entrez) expression.")

        dataSet$data.anot <- RemoveDuplicates(data.anot, lvlOpt, quiet = F)
        dataSet$id.current <- "entrez"
        dataSet$annotated <- T
      } else {
        current.msg <<- paste(current.msg, "No gene level summarization was performed.")
      }
    }
  } else { # no conversion will be performed
    feature.vec <- rownames(data.proc)
    anot.id <- feature.vec
    perct <- 100
    hit.inx <- !is.na(anot.id)
    matched.len <- length(feature.vec) # dummies
    1
    current.msg <<- paste("No annotation was performed. Make sure organism and gene ID are specified correctly!")
  }
  # need to save the ids (mixed gene annotation and original id)
  # in case, users needs to keep unannotated features
  # this need to be updated to gether with data from now on
  dataSet$data.norm <- dataSet$data.anot
  dataSet <<- dataSet

  saveRDS(dataSet$data.anot, file = "orig.data.anot") # keep original copy, not in mem

  totalCount <- sum(colSums(dataSet$data.anot))
  avgCount <- sum(colSums(dataSet$data.anot)) / ncol(dataSet$data.anot)
  minCount <- min(colSums(dataSet$data.anot))
  maxCount <- max(colSums(dataSet$data.anot))

  if (length(dataSet$meta.info) == 1) {
    lvls <- paste(levels(dataSet$meta.info[, 1]), collapse = "; ")
  } else {
    conc1 <- paste0("<b>", colnames(dataSet$meta.info)[1], "</b>", ": ", paste(levels(dataSet$meta.info[, 1]), collapse = "; "))
    conc2 <- paste0("<b>", colnames(dataSet$meta.info)[2], "</b>", ": ", paste(levels(dataSet$meta.info[, 2]), collapse = "; "))
    lvls <- paste("Two factors found -", conc1, conc2)
  }
  summaryVec <<- c(matched.len, perct, length(anot.id), sum(!hit.inx), ncol(dataSet$data.anot), ncol(dataSet$meta.info), sprintf("%4.2e", signif(totalCount, 3)), sprintf("%4.2e", signif(avgCount, 3)), sprintf("%4.2e", signif(minCount, 3)), sprintf("%4.2e", signif(maxCount, 3)), lvls)
  return(matched.len)
}


PerformExpressNormalization <- function(norm.opt, var.thresh, count.thresh, abundance, filterUnmapped) {
  print("normalizing ....")
  msg <- "Only features with annotations are kept for further analysis."

  if (filterUnmapped == "false") {
    # need to update those with annotations
    data1 <- readRDS("data.proc.rds")
    anot.id <- readRDS("annotation.rds")
    hit.inx <- !is.na(anot.id)
    rownames(data1)[hit.inx] <- anot.id[hit.inx]
    data1 <- RemoveDuplicates(data1, "mean", quiet = T)
    raw.data.anot <- data <- dataSet$data.anot <- data1
  } else {
    raw.data.anot <- data <- readRDS("orig.data.anot")
  }

  if (dataSet$type == "count") {
    sum.counts <- apply(data, 1, sum, na.rm = TRUE)
    rm.inx <- sum.counts < count.thresh
    data <- data[!rm.inx, ]
    msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low counts.", collapse = " ")
  } else {
    avg.signal <- apply(data, 1, mean, na.rm = TRUE)
    p05 <- quantile(avg.signal, 0.05)
    nrow(data)
    rm.inx <- avg.signal < p05
    data <- data[!rm.inx, ]
    msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low relative abundance (average expression signal).", collapse = " ")
  }

  data <- PerformDataNormalization(data, norm.opt)
  if (length(data) == 1 && data == 0) {
    return(0)
  }
  msg <- paste(norm.msg, msg)

  filter.val <- apply(data, 1, IQR, na.rm = T)
  "Interquantile Range"
  rk <- rank(-filter.val, ties.method = "random")
  kp.pct <- (100 - var.thresh) / 100

  remain <- rk < nrow(data) * kp.pct
  data <- data[remain, ]
  msg <- paste(msg, paste("Filtered ", sum(!remain), " low variance genes based on IQR"), collapse = " ")

  dataSet$data.anot <- raw.data.anot[remain, ]
  dataSet$data.norm <- data

  # save normalized data for download user option
  write.csv(dataSet$data.norm, file = "data_normalized.csv")

  current.msg <<- msg
  dataSet <<- dataSet

  return(1)
}


# note, we do both filtering and normalization
PerformDataNormalization <- function(data, norm.opt) {
  set.seed(1337)
  msg <- NULL
  row.nms <- rownames(data)
  col.nms <- colnames(data)
  if (norm.opt == "log") {
    min.val <- min(data[data > 0], na.rm = T) / 10
    numberOfNeg <- sum(data <= 0, na.rm = TRUE) + 1
    totalNumber <- length(data)
    if ((numberOfNeg / totalNumber) > 0.2) {
      msg <- paste(msg, "Can't perform log2 normalization, over 20% of data are negative. Try a different method or maybe the data already normalized?", collapse = " ")
      print(msg)
      norm.msg <<- current.msg <<- msg
      return(0)
    }
    data[data <= 0] <- min.val
    data <- log2(data)
    msg <- paste(msg, "Log2 transformation.", collapse = " ")
  } else if (norm.opt == "vsn") {
    require(limma)
    data <- normalizeVSN(data)
    msg <- paste(msg, "VSN normalization.", collapse = " ")
  } else if (norm.opt == "quantile") {
    require("preprocessCore")
    data <- normalize.quantiles(data, copy = TRUE)
    msg <- paste(msg, "Quantile normalization.", collapse = " ")
  } else if (norm.opt == "combined") {
    require(limma)
    data <- normalizeVSN(data)
    require("preprocessCore")
    data <- normalize.quantiles(data, copy = TRUE)
    msg <- paste(msg, "VSN followed by quantile normalization.", collapse = " ")
  } else if (norm.opt == "logcount") { # for count data, do it in DE analysis, as it is dependent on design matrix
    require(edgeR)
    nf <- calcNormFactors(data)
    y <- voom(data, plot = F, lib.size = colSums(data) * nf)
    data <- y$E # copy per million
    msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse = " ")
  } else {
    # should do best guess for count data for plotting and filtering
    if (dataSet$type == "count") {
      if (sum(data > 100) > 100) { # now we think it is raw counts
        require(edgeR)
        nf <- calcNormFactors(data)
        y <- voom(data, plot = F, lib.size = colSums(data) * nf)
        data <- y$E # copy per million
      }
    }
    msg <- paste(msg, "No log normalization was performed.", collapse = " ")
    print(msg)
  }
  norm.msg <<- msg
  rownames(data) <- row.nms
  colnames(data) <- col.nms
  return(data)
}

# note, setup the main class, keep the original order
SetMainClass <- function(cls.lbl) {
  lbls <- as.character(dataSet$meta.info[[cls.lbl]])
  lvls.orig <- unique(lbls)
  cls <- factor(lbls, levels = lvls.orig, ordered = T)
  dataSet$cls <- cls # record main cls
  dataSet <<- dataSet
  return(levels(cls))
}

SetSelectedMetaInfo <- function(meta0, meta1, block1) {
  if (meta0 == "NA") {
    return(0)
  } else {
    cls <- dataSet$meta.info[, meta0]
    dataSet$fst.cls <- cls # for PCA plotting
    block <- NULL
    dataSet$sec.cls <- "NA"
    if (meta1 != "NA") {
      if (block1) {
        block <- dataSet$meta.info[, meta1]
      } else { # two factor
        cls <- interaction(dataSet$meta.info[, c(meta0, meta1)], sep = ".", lex.order = TRUE)
      }
      dataSet$sec.cls <- dataSet$meta.info[, meta1] # for pca coloring
    }
    dataSet$cls <- cls # record main cls;
    dataSet$block <- block
    dataSet <<- dataSet
    return(levels(cls))
  }
}

SetupDesignMatrix <- function(deMethod) {
  cls <- dataSet$cls
  design <- model.matrix(~ 0 + cls) # no intercept
  colnames(design) <- levels(cls)
  dataSet$design <- design
  dataSet$de.method <- deMethod
  dataSet <<- dataSet

  return(1)
}

# perform differential analysis
# default: all pair-wise comparison (A-B) + (B-C) + (A-C)
# custom: only compare two groups (A-C)
# time: compare consecutive groups (B-A) + (C-B)
# reference: all others against common reference (A-C) + (B-C)
# nested: (A-B)+(C-D)
PerformDEAnal <- function(anal.type = "default", par1 = NULL, par2 = NULL, nested.opt = "intonly") {
  set.seed(1337)
  par1 <<- par1
  par2 <<- par2
  nested.opt <<- nested.opt
  myargs <- list()
  cls <- dataSet$cls
  dataSet$comp.type <- anal.type
  grp.nms <- levels(cls)

  if (anal.type == "default") {
    inx <- 0
    for (m in 1:(length(grp.nms) - 1)) {
      for (n in (m + 1):length(grp.nms)) {
        inx <- inx + 1
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep = "")
      }
    }
    filename <- "SigGene_pairwise"
  } else if (anal.type == "time") {
    for (i in 2:length(grp.nms)) {
      myargs[[i - 1]] <- paste(grp.nms[i], "-", grp.nms[i - 1], sep = "")
    }
    filename <- "SigGene_time_series"
  } else if (anal.type == "custom") {
    grp.nms <- strsplit(par1, " vs. ")[[1]]
    myargs[[1]] <- paste(grp.nms, collapse = "-")
    dataSet$grp.nms <- grp.nms
    filename <- paste("SigGene_", paste(grp.nms, collapse = "_vs_"), sep = "")
  } else if (anal.type == "reference") {
    ref <- par1
    cntr.cls <- grp.nms[grp.nms != ref]
    myargs <- as.list(paste(cntr.cls, "-", ref, sep = ""))
    filename <- paste("SigGene_reference_", ref, sep = "")
  } else if (anal.type == "nested") {
    grp.nms1 <- strsplit(par1, " vs. ")[[1]]
    grp.nms2 <- strsplit(par2, " vs. ")[[1]]
    if (all(grp.nms1 == grp.nms2)) {
      current.msg <<- paste("The two nested groups are the same. Please choose two different groups.")
      return(0)
    }
    grp.nms <- unique(c(grp.nms1, grp.nms2))
    if (nested.opt == "intonly") {
      myargs[[1]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    } else {
      myargs[[1]] <- paste(grp.nms1, collapse = "-")
      myargs[[2]] <- paste(grp.nms2, collapse = "-")
      myargs[[3]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    }
    filename <- paste("SigGene_nested_", paste(paste(grp.nms1, collapse = "_vs_"), "_", paste(grp.nms2, collapse = "_vs_"), sep = ""), sep = "")
  } else {
    print(paste("Not supported: ", anal.type))
  }

  library(limma)
  design <- dataSet$design
  myargs[["levels"]] <- design
  contrast.matrix <- do.call(makeContrasts, myargs)
  if (dataSet$de.method == "limma") {
    if (is.null(dataSet$block)) {
      fit <- lmFit(dataSet$data.norm, design)
    } else {
      corfit <- duplicateCorrelation(dataSet$data.norm, design, block = dataSet$block)
      fit <- lmFit(dataSet$data.norm, design, block = dataSet$block, correlation = corfit$consensus)
    }

    if (!is.fullrank(design)) {
      current.msg <<- paste("This metadata combination is not full rank! Please use other combination.")
      return(0)
    }

    df.residual <- fit$df.residual
    if (all(df.residual == 0)) {
      current.msg <<- paste("There is not enough replicates in each group (no residual degrees of freedom)!")
      return(0)
    }
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    topFeatures <- topTable(fit2, number = Inf, adjust.method = "fdr")
  } else if (dataSet$de.method == "deseq2") {
    # only for small data set (< 60)
    if (length(cls) > 60) {
      current.msg <<- "For large sample size (>60), use limma or edgeR."
      return(0)
    } else { # use microservice
      deseq.in <- list(dataSet = dataSet, contrast.matrix = contrast.matrix)
      saveRDS(deseq.in, "deseq_in.rds")
      return(1)
    }
  } else {
    library(edgeR)
    y <- DGEList(counts = dataSet$data.anot, group = dataSet$cls)
    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y, design, verbose = FALSE)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, contrast = contrast.matrix)
    topFeatures <- topTags(lrt, n = Inf)$table
    nms <- colnames(topFeatures)
    nms[which(nms == "FDR")] <- "adj.P.Val"
    colnames(topFeatures) <- nms
  }
  dataSet$filename <- filename
  dataSet$resTable <- topFeatures
  dataSet <<- dataSet
  return(1)
}

# update result based on new cutoff
GetSigGenes <- function(res.nm, p.lvl, fc.lvl, update = T, inx) {
  total <- nrow(dataSet$resTable)
  resTable <- dataSet$resTable
  filename <- dataSet$filename
  filename <- paste(filename, "_", res.nm, ".csv", sep = "")
  if (update) {
    current.msg <<- ""
  }
  # select based on p-value
  if (dataSet$type == "array") {
    hit.inx.p <- resTable$adj.P.Val <= p.lvl
  } else {
    hit.inx.p <- resTable$adj.P.Val <= p.lvl
  }

  resTable <- resTable[hit.inx.p, , drop = F]
  if (nrow(resTable) == 0) {
    current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff.")
  }
  # now rank by logFC, note, the logFC for each comparisons
  # are returned in resTable before the AveExpr columns
  # for two-class, only one column, multiple columns can be involved
  # for > comparisons - in this case, use the largest logFC among all comparisons
  # if (fc.lvl > 0){ # further filter by logFC
  if (dataSet$de.method == "limma") {
    hit.inx <- which(colnames(resTable) == "AveExpr")
  } else if (dataSet$de.method == "deseq2") {
    hit.inx <- which(colnames(resTable) == "baseMean")
  } else {
    hit.inx <- which(colnames(resTable) == "logCPM")
  }
  maxFC.inx <- hit.inx - 1 # not sure if this is also true for edgeR
  logfc.mat <- resTable[, 1:maxFC.inx, drop = F]
  pos.mat <- abs(logfc.mat)
  fc.vec <- apply(pos.mat, 1, max)
  hit.inx.fc <- fc.vec >= fc.lvl
  resTable <- resTable[hit.inx.fc, , drop = F]
  if (nrow(resTable) == 0) {
    current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff.")
  }
  # }

  ### Note, rowname of resTable must be entrez ID

  de.Num <- nrow(resTable)

  # display at most 5000 genes for the server (two main reasons)
  # 1) should not have more 22% (human: 23000) DE of all genes (biological)
  # 2) IE canvas can display no more than 6800 pixels (computational)
  if (nrow(resTable) > 5000) {
    resTable <- resTable[1:5000, ]
    current.msg <<- paste(current.msg, " Due to computational constraints, only top 5000 genes will be used. ", collapse = "\n")
  }

  # may need to update data, class and meta.info
  data <- dataSet$data.norm
  cls <- dataSet$cls
  meta.info <- dataSet$meta.info
  grp.nms <- levels(cls)

  hit.inx <- cls %in% grp.nms
  if (sum(hit.inx) < length(hit.inx)) {
    current.msg <<- paste(current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse = ", "), "are included.")
    cls <- factor(cls[hit.inx])
    levels(cls)
    data <- data[, hit.inx]
    meta.info <- dataSet$meta.info[hit.inx, ]
  }

  saveRDS(data, file = "data.stat")
  dataSet$resTable <- dataSet$resTable[order(dataSet$resTable$adj.P.Val), ]
  dataSet$resTable <- dataSet$resTable[which(!rownames(dataSet$resTable) %in% rownames(resTable)), ]
  dataSet$resTable <- rbind(resTable, dataSet$resTable)

  dataSet$sig.mat <- resTable
  if (dataSet$annotated) { # annotated to entrez
    anot.id <- rownames(dataSet$resTable)
    gene.anot <- doEntrezIDAnot(anot.id)
    write.csv(cbind(EntrezID = anot.id, signif(dataSet$resTable, 5), Symbols = gene.anot$symbol, Name = gene.anot$name), row.names = F, file = filename)
  } else if (file.exists("annotation.rds")) { # annotation information available
    anot.id <- readRDS("annotation.rds")
    feature.vec <- rownames(dataSet$resTable)
    entrez.vec <- anot.id[feature.vec]
    gene.anot <- doEntrezIDAnot(entrez.vec)
    write.csv(cbind(signif(dataSet$resTable, 5), EntrezID = entrez.vec, Symbols = gene.anot$symbol, Name = gene.anot$name), row.names = F, file = filename)
    rownames(gene.anot) <- feature.vec
  } else {
    gene.anot <- NULL
    write.csv(signif(resTable, 5), file = filename)
  }
  if (is.null(gene.anot)) {
    dataSet$sig.genes.symbols <- rep("NA", nrow(resTable))
  } else {
    dataSet$sig.genes.symbols <- gene.anot$symbol
  }
  dataSet$cls.stat <- cls
  dataSet$meta.stat <- meta.info

  # now do protein mapping for network only applicable for annotated

  dataSet$name <- res.nm

  gene <- rownames(resTable)

  logFC <- unname(logfc.mat[, 1])
  geneList <- paste(gene, logFC, collapse = "\n")
  up <- nrow(resTable[which(logfc.mat[, selectedFactorInx] > fc.lvl), ])
  down <- nrow(resTable[which(logfc.mat[, selectedFactorInx] < -fc.lvl), ])

  dataSet <<- dataSet
  data.norm <- dataSet$data.norm
  colnames(data.norm) <- NULL
  lst <- list(colnames(dataSet$data.norm), data.norm, dataSet$meta.info, dataSet$resTable, rownames(data.norm), org = data.org)
  library(RJSONIO)
  json.obj <- toJSON(lst)
  sink("NetworkAnalyst_matrix.json")
  cat(json.obj)
  return(c(filename, de.Num, geneList, total, up, down))
}

GetExpressResultColNames <- function() {
  resT <- readRDS("ExpressResT.rda")
  colnames(resT)
}

GetExpressResultGeneIDs <- function() {
  return(rownames(dataSet$resTable))
}

GetExpressGeneIDType <- function() {
  return(dataSet$id.current)
}

GetExpressResultMatrix <- function(inxt) {
  inxt <- as.numeric(inxt)
  if (dataSet$de.method == "limma") {
    inx <- match("AveExpr", colnames(dataSet$resTable))
  } else if (dataSet$de.method == "deseq2") {
    inx <- match("baseMean", colnames(dataSet$resTable))
  } else {
    inx <- match("logCPM", colnames(dataSet$resTable))
  }
  res <- dataSet$resTable
  res <- res[, -(1:inx - 1)]
  res <- cbind(dataSet$resTable[, inxt], res)
  colnames(res)[1] <- colnames(dataSet$resTable)[inxt]
  saveRDS(res, "ExpressResT.rda")
  return(signif(as.matrix(res), 5))
}


GetExpressResultGeneIDLinks <- function() {
  ids <- rownames(dataSet$resTable)
  symbs <- doEntrez2SymbolMapping(ids)
  annots <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids, "' target='_blank'>", symbs, "</a>", sep = "")
  return(annots)
}

GetExpressResultGeneSymbols <- function() {
  return(dataSet$sig.genes.symbols)
}

GetFactorNb <- function() {
  return(length(dataSet$meta.info))
}

PlotDataBox <- function(boxplotName, dpi, format) {
  qc.boxplot(dataSet$data.norm, boxplotName, dpi, format)
}

PlotDataPCA <- function(pcaName, dpi, format, factor) {
  qc.pcaplot(dataSet$data.norm, pcaName, dpi, format, factor)
}

PlotDataMeanStd <- function(densityName, dpi, format) {
  qc.meanstd(dataSet$data.norm, densityName, dpi, format)
}


qc.density <- function(imgNm, dpi = 72, format, factor) {
  library("ggplot2")
  dataSet$data.norm
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  dpi <- as.numeric(dpi)

  df <- data.frame(dataSet$data.norm, stringsAsFactors = FALSE)
  df <- stack(df)
  sampleNms <- gsub("-", ".", colnames(dataSet$data.norm))
  if (length(dataSet$meta.info) == 2) {
    Factor1 <- as.vector(dataSet$meta.info[, 1])
    factorNm1 <- colnames(dataSet$meta.info)[1]
    conv <- data.frame(ind = sampleNms, Factor1 = Factor1)
    colnames(conv) <- c("ind", factorNm1)
    df1 <- merge(df, conv, by = "ind")
    Factor2 <- as.vector(dataSet$meta.info[, 2])
    factorNm2 <- colnames(dataSet$meta.info)[2]
    conv <- data.frame(ind = sampleNms, Factor2 = Factor2)
    colnames(conv) <- c("ind", factorNm2)
    df1 <- merge(df1, conv, by = "ind")
    df2 <- melt(df1, measure.vars = c(factorNm1, factorNm2))
    colnames(df2)[4] <- "Conditions"
    g <- ggplot(df2, aes(x = values)) + geom_line(aes(color = Conditions, group = ind), stat = "density", alpha = 0.6) + facet_grid(. ~ variable)
    width <- 12
    height <- 6
  } else {
    Conditions <- as.character(dataSet$meta.info[, 1])
    conv <- data.frame(ind = sampleNms, Conditions = Conditions)
    df1 <- merge(df, conv, by = "ind")
    g <- ggplot(df1, aes(x = values)) + geom_line(aes(color = Conditions, group = ind), stat = "density", alpha = 0.6)
    width <- 8
    height <- 6
  }
  Cairo(file = imgNm, width = width, height = height, type = format, bg = "white", dpi = dpi, unit = "in")
  print(g)
  dev.off()
}

qc.densitySample <- function(dat, imgNm, dpi = 72, format) {
  library("ggplot2")

  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  dpi <- as.numeric(dpi)
  Cairo(file = imgNm, width = 8, height = 6, type = format, bg = "white", dpi = dpi, unit = "in")
  df <- data.frame(dataSet$data.norm, stringsAsFactors = FALSE)
  df <- stack(df)
  g <- ggplot(df, aes(x = values, color = ind)) +
    geom_density()

  print(g)
  dev.off()
}



qc.meanstd <- function(dat, imgNm, dpi = 72, format = "png") {
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  Cairo(file = imgNm, width = 8, height = 6, type = format, bg = "white", dpi = dpi, unit = "in")
  meanSdPlot(dat, ranks = FALSE)
  dev.off()
}

qc.boxplot <- function(dat, imgNm, dpi = 72, format = "png") {
  dpi <- as.numeric(dpi)
  library("ggplot2")
  library("lattice")
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  subgene <- 10000
  if (nrow(dat) > subgene) {
    set.seed(28051968)
    sg <- sample(nrow(dat), subgene)
    Mss <- dat[sg, , drop = FALSE]
  } else {
    Mss <- dat
  }

  subsmpl <- 100
  if (ncol(Mss) > subsmpl) {
    set.seed(28051968)
    ss <- sample(ncol(Mss), subsmpl)
    Mss <- Mss[, ss, drop = FALSE]
  } else {
    Mss <- Mss
  }

  sample_id <- rep(seq_len(ncol(Mss)), each = nrow(Mss))
  values <- as.numeric(Mss)

  df <- cbind(values, sample_id)

  df <- data.frame(df)
  df$sample_id <- factor(df$sample_id)
  xlower <- unname(quantile(df$values, probs = c(0.01, 0.99), na.rm = TRUE)[1])
  xupper <- unname(quantile(df$values, probs = c(0.01, 0.99), na.rm = TRUE)[2])
  height <- length(unique(df$sample_id)) * 20
  if (height < 450) {
    height <- 450
  }
  bp <- ggplot(df, aes(sample_id, values)) +
    ylab("Values") + xlab("Samples") + scale_x_discrete(labels = colnames(dataSet$data.norm)) + ylim(xlower, xupper) + stat_boxplot(geom = "errorbar", color = "black") + geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4)
  bp <- bp + coord_flip()

  Cairo(file = imgNm, width = 600 * dpi / 72, height = height * dpi / 72, unit = "px", dpi = dpi, type = format, bg = "white")
  print(bp)
  dev.off()
}

qc.pcaplot <- function(x, imgNm, dpi = 72, format = "png", factor) {
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  library("lattice")
  library("ggplot2")
  library("reshape")
  pca <- prcomp(t(na.omit(x)))
  imp.pca <- summary(pca)$importance
  xlabel <- paste0("PC1", " (", 100 * round(imp.pca[2, ][1], 3), "%)")
  ylabel <- paste0("PC2", " (", 100 * round(imp.pca[2, ][2], 3), "%)")
  names <- colnames(x)
  pca.res <- as.data.frame(pca$x)
  pca.res <- pca.res[, c(1, 2)]
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1)
  ylim <- GetExtendRange(pca.res$PC2)
  if (length(dataSet$meta.info) == 2) {
    Factor1 <- as.vector(dataSet$meta.info[, 1])
    factorNm1 <- colnames(dataSet$meta.info)[1]
    pca.res[, factorNm1] <- Factor1
    Factor2 <- as.vector(dataSet$meta.info[, 2])
    factorNm2 <- colnames(dataSet$meta.info)[2]
    pca.res[, factorNm2] <- Factor2
    pca.rest <- melt(pca.res, measure.vars = c(factorNm1, factorNm2))
    colnames(pca.rest)[4] <- "Conditions"
    pca.rest$names <- c(rownames(pca.res), rownames(pca.res))
    if (length(pca.rest$names) > 20) {
      pcafig <- ggplot(pca.rest, aes(x = PC1, y = PC2, color = Conditions, label = pca.rest$names)) +
        geom_point(size = 3, alpha = 0.5) + xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) + facet_grid(. ~ variable)
    } else {
      require("ggrepel")
      pcafig <- ggplot(pca.rest, aes(x = PC1, y = PC2, color = Conditions, label = pca.rest$names)) +
        geom_point(size = 4) + xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) + geom_text_repel(force = 1.5) + facet_grid(. ~ variable)
    }
    width <- 12
    height <- 6
  } else {
    Factor <- dataSet$meta.info[, 1]
    pca.rest <- pca.res
    pca.rest$Conditions <- Factor
    pca.rest$names <- rownames(pca.res)
    if (length(rownames(pca.res)) > 20) {
      pcafig <- ggplot(pca.rest, aes(x = PC1, y = PC2, color = Conditions)) +
        geom_point(size = 3, alpha = 0.5) + xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel)
    } else {
      require("ggrepel")
      pcafig <- ggplot(pca.rest, aes(x = PC1, y = PC2, color = Conditions, label = rownames(pca.res))) +
        geom_point(size = 4) + xlim(xlim) + ylim(ylim) + xlab(xlabel) + ylab(ylabel) + geom_text_repel(force = 1.5)
    }
    width <- 8
    height <- 6
  }

  Cairo(file = imgNm, width = width, height = height, type = format, bg = "white", unit = "in", dpi = dpi)
  print(pcafig)
  dev.off()
}


PlotLibSizeView <- function(imgNm, dpi = 72, format = "png", factor) {
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  dpi <- as.numeric(dpi)
  data_bef <- data.matrix(dataSet$data.anot)

  smpl.sums <- colSums(data_bef)

  library("ggplot2")
  data_bef <- data.matrix(dataSet$data.anot)
  smpl.sums <- colSums(data_bef)
  names(smpl.sums) <- colnames(data_bef)
  sampleNms <- names(smpl.sums)
  df <- data.frame(count = smpl.sums, ind = colnames(data_bef))

  if (length(dataSet$meta.info) == 2) {
    Factor1 <- as.vector(dataSet$meta.info[, 1])
    factor1Nm <- colnames(dataSet$meta.info)[1]
    conv <- data.frame(ind = sampleNms, Factor1 = Factor1)
    colnames(conv) <- c("ind", factor1Nm)
    df1 <- merge(df, conv, by = "ind")
    Factor2 <- as.vector(dataSet$meta.info[, 2])
    factor2Nm <- colnames(dataSet$meta.info)[2]
    conv <- data.frame(ind = sampleNms, Factor2 = Factor2)
    colnames(conv) <- c("ind", factor2Nm)
    df1 <- merge(df1, conv, by = "ind")
    df2 <- melt(df1, measure.vars = c(factor1Nm, factor2Nm))
    colnames(df2)[4] <- "Conditions"
    if (length(df2$ind) > 20) {
      g <- ggplot(df2, aes(x = Conditions, y = count, fill = Conditions)) +
        geom_dotplot(
          binaxis = "y", stackdir = "center",
          position = position_dodge(), dotsize = 0.7
        ) + ylab("Sum") + facet_grid(. ~ variable)
    } else {
      g <- ggplot(df2, aes(x = Conditions, y = count, fill = Conditions, label = ind)) +
        geom_dotplot(
          binaxis = "y", stackdir = "center",
          position = position_dodge(), dotsize = 0.7
        ) + geom_text_repel(force = 5) + ylab("Sum") + facet_grid(. ~ variable)
    }
    width <- 12
    height <- 6
  } else {
    Conditions <- as.character(dataSet$meta.info[, 1])
    conv <- data.frame(ind = sampleNms, Conditions = Conditions)
    df1 <- merge(df, conv, by = "ind")
    if (length(df1$ind) > 20) {
      g <- ggplot(df1, aes(x = Conditions, y = count, fill = Conditions)) +
        geom_dotplot(
          binaxis = "y", stackdir = "center",
          position = position_dodge(), dotsize = 0.7
        ) + xlab("Sum")
    } else {
      g <- ggplot(df1, aes(x = Conditions, y = count, label = ind, fill = Conditions)) +
        geom_dotplot(
          binaxis = "y", stackdir = "center",
          position = position_dodge(), dotsize = 0.7
        ) + geom_text_repel(force = 5) + xlab("Sum")
    }
    width <- 8
    height <- 6
  }

  Cairo(file = imgNm, width = width, height = height, unit = "in", type = format, bg = "white", dpi = dpi)

  # names(smpl.sums) <- colnames(data_bef);
  print(g)
  dev.off()
}

PlotMDS <- function(imgName, format) {
  library(edgeR)
  library(RColorBrewer)
  imgName <- paste(imgName, "dpi", dpi, ".", format, sep = "")
  Cairo(file = imgName, width = 580, type = format, bg = "white", dpi = 72)
  levels(dataSet$cls) <- brewer.pal(nlevels(dataSet$cls), "Set1")
  col.group <- dataSet$cls
  col.group <- as.character(col.group)
  plotMDS(dataSet$data.norm, col = col.group, xlab = "Dimension 2", ylab = "Dimension 1")
  title(main = "MDS")
  dev.off()
}

GetSummaryData <- function() {
  return(summaryVec)
}

GetDensityPlot <- function() {
  data <- list()
  dat <- dataSet$data.norm
  for (i in 1:ncol(dat)) {
    data[[i]] <- dat[, i]
    names(data)[i] <- colnames(dat)[i]
  }
  densityList <- list()
  for (i in 1:length(data)) {
    d <- density(data[[i]])
    df <- data.frame(d$x, d$y)
    colnames(df) <- c("x", "y")
    densityList[[i]] <- df
  }
  names(densityList) <- colnames(dataSet$data.norm)
  jsonNm <- "density.json"
  lst <- list()

  lst <- list(
    density = densityList,
    class = dataSet$cls
  )
  library(RJSONIO)
  json.obj <- toJSON(lst)
  sink(jsonNm)
  cat(json.obj)
}

PlotMAPlot <- function(imgNm, dpi = 72, format, pvalue, fc, inx) {
  library("ggplot2")
  dpi <- as.numeric(dpi)
  inx <- as.numeric(inx)
  pvalue <- as.numeric(pvalue)
  fc <- as.numeric(fc)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  Cairo(file = imgNm, width = 700, height = 560, type = format, bg = "white", dpi = 100)
  pvalue <- as.numeric(pvalue)
  fc <- as.numeric(fc)
  res <- dataSet$resTable
  if (dataSet$de.method == "deseq2") {
    # res = res[!apply(sapply(res, function(x) abs(scale(x)) >= 5), 1, any), ];
    res["log2(baseMean)"] <- log2(res$baseMean)
  }

  # select based on p-value
  if (dataSet$type == "array") {
    res$significant <- ifelse(res[, inx] > fc & res$P.Value < pvalue, 1, ifelse(res[, inx] < -fc & res$P.Value < pvalue, -1, 0))
  } else {
    res$significant <- ifelse(res[, inx] > fc & res$adj.P.Val < pvalue, 1, ifelse(res[, inx] < -fc & res$adj.P.Val < pvalue, -1, 0))
  }
  res$significant <- as.factor(res$significant)
  yCol <- colnames(res)[inx]

  if (dataSet$de.method == "limma") {
    maplot <- ggplot(res, aes_string(x = "AveExpr", y = yCol, color = "significant"))
  } else if (dataSet$de.method == "deseq2") {
    maplot <- ggplot(res, aes_string(x = "log2(baseMean)", y = yCol, color = "significant"))
  } else {
    maplot <- ggplot(res, aes_string(x = "logCPM", y = yCol, color = "significant"))
  }

  maplot <- maplot +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(color = "blue3", yintercept = 0) +
    stat_smooth(se = FALSE, method = "loess", color = "red3") +
    scale_color_manual(values = c("-1" = "green", "0" = "black", "1" = "red")) + theme(legend.position = "none")
  print(maplot)
  dev.off()
}

GetMetaColLength <- function() {
  if (dataSet$de.method == "limma") {
    inx <- match("AveExpr", colnames(dataSet$resTable))
  } else if (dataSet$de.method == "deseq2") {
    inx <- match("baseMean", colnames(dataSet$resTable))
  } else {
    inx <- match("logCPM", colnames(dataSet$resTable))
  }
  resT <- dataSet$resTable
  resT <- resT[, 1:inx - 1]
  return(length(colnames(resT)))
}


GetMetaInfoLength <- function() {
  return(length(dataSet$meta.info))
}

GetMetaCol <- function() {
  colNms <- colnames(dataSet$resTable)
  if (dataSet$de.method == "limma") {
    inx <- match("AveExpr", colNms)
  } else if (dataSet$de.method == "deseq2") {
    inx <- match("baseMean", colNms)
  } else {
    inx <- match("logCPM", colNms)
  }
  resT <- dataSet$resTable
  resT <- resT[, 1:inx - 1]
  return(colnames(resT))
}

Volcano.Anal <- function(paired = FALSE, fcthresh, threshp, analType, inx) {
  inx <- as.numeric(inx)
  print("Prepare volcano anal")
  if (anal.type == "metadata") {
    if (dataSet$name != selDataNm) {
      dataSet <- readRDS(selDataNm)
    }
    data <- as.matrix(inmex.ind[selDataNm][[1]])
    p.value <- data[, "Pval"]
    fcthresh <- 0
  } else {
    data <- as.matrix(dataSet$resTable)

    if (dataSet$type == "array") {
      p.value <- data[, "adj.P.Val"]
    } else {
      p.value <- data[, "adj.P.Val"]
    }
  }
  fcthreshu <<- fcthresh

  if (analType == "qPCR") {
    inx.p <- p.value < 1
  } else {
    inx.p <- p.value <= threshp
  }
  zero.inx <- p.value == 0
  if (sum(zero.inx) > 0) {
    p.value[zero.inx] <- min(p.value[!zero.inx]) / 10
  }
  p.log <- -log10(p.value)

  if (dataSet$annotated) { # annotated to entrez
    anot.id <- rownames(data)
    gene.anot <- doEntrezIDAnot(anot.id)
  } else {
    anot.id <- rownames(data)
    gene.anot <- data.frame(gene_id = anot.id, symbol = anot.id, stringsAsFactors = FALSE)
    init.lib <<- "NA"
  }

  # gene symbol to be used for boxplot

  # create a named matrix of sig vars for display
  fc.log <- data[, inx]
  hit.maxPos <- (which(fc.log > 10))
  hit.maxNeg <- (which(fc.log < -10))
  fc.log[hit.maxPos] <- 10
  fc.log[hit.maxNeg] <- 10
  # fc.all <- res$fc.all;

  cs_1 <- p.value < threshp
  if (fcthresh != 0) {
    inx.up <- fc.log > fcthresh & cs_1
    inx.down <- fc.log < -fcthresh & cs_1
  } else {
    inx.up <- fc.log > 0 & cs_1
    inx.down <- fc.log < 0 & cs_1
  }

  # create named sig table for display
  inx.imp <- (inx.up | inx.down) & inx.p
  sig.var <- cbind(fc.log[inx.imp, drop = F], p.value[inx.imp, drop = F], p.log[inx.imp, drop = F])
  colnames(sig.var) <- c("log2(FC)", "p.value", "-log10(p)")
  # first order by log(p), then by log(FC)
  ord.inx <- order(sig.var[, 3], abs(sig.var[, 1]), decreasing = T)
  sig.var <- sig.var[ord.inx, , drop = F]

  sig.var <- signif(sig.var, 5)
  sig.var1 <- sig.var
  sig.var1 <- cbind(rownames(sig.var), sig.var)
  colnames(sig.var1) <- c("name", "log2(FC)", "p.value", "-log10(p)")

  ###########################
  ## for Volcano data
  ##########################

  if (init.lib != "NA") {
    PerformVolcanoEnrichment("abc", init.lib, "null", "all", inx)
  }

  fileName <- "volcano.csv"
  jsonNm <- "volcano.json"
  library(RJSONIO)
  json.obj <- toJSON(sig.var1)
  sink(jsonNm)
  cat(json.obj)
  sink()
  write.csv(signif(sig.var, 5), file = fileName)
  colnames(gene.anot)[1] <- "anot.id"
  volcano <- list(
    raw.threshx = fcthresh,
    raw.threshy = threshp,
    paired = paired,
    thresh.y = -log10(threshp),
    fc.symb = rownames(data),
    fc.log = fc.log,
    fc.log.uniq = jitter(fc.log),
    inx.up = inx.up,
    inx.down = inx.down,
    p.log = p.log,
    inx.p = inx.p,
    sig.mat = sig.var,
    conv = gene.anot
  )

  library(RJSONIO)
  json.obj <- toJSON(volcano)
  sink("volcano2.json")
  cat(json.obj)
  sink()

  if (init.lib == "NA") {
    enr.mat <- "NA"
  }
  write.csv(enr.mat, file = "enrichment_result.csv", row.names = T)
  sink("enrichment_result.json")
  cat(json.obj)
  sink()
}

# perform limma on given two groups selected
# used by integarative analysis
PerformLimmaDE <- function(dataName, grps, p.lvl, fc.lvl = NULL) {
  print("doing differential analysis ....")
  dataSet <- readRDS(dataName)
  dataSet$pval <- p.lvl
  if (length(levels(dataSet$cls)) > 2) {
    grp.nms <- strsplit(grps, " vs. ")[[1]]
    sel.inx <- as.character(dataSet$cls) %in% grp.nms
  } else {
    sel.inx <- rep(T, ncol(dataSet$data))
  }

  group <- factor(dataSet$cls[sel.inx]) # note regenerate factor to drop levels
  data <- dataSet$data[, sel.inx]

  res.limma <- PerformLimma(data, group)
  res.all <- GetLimmaResTable(res.limma$fit.obj)

  if (!is.null(fc.lvl)) {
    hit.inx <- abs(res.all$logFC) >= fc.lvl & res.all$adj.P.Val <= p.lvl
  } else {
    hit.inx <- res.all$adj.P.Val <= p.lvl
  }
  # note, hit.inx can contain NA, not T/F
  hit.inx <- which(hit.inx)
  res <- res.all[hit.inx, ]

  # rm .txt suffix for new names
  shortNm <- substring(dataName, 0, nchar(dataName) - 4)
  write.csv(signif(res[, -1], 5), file = paste("SigGenes_", shortNm, ".csv", sep = ""))

  sig.count <- nrow(res)
  rownames(res)
  res.mat <- cbind(res.all$logFC, res.all$adj.P.Val)
  rownames(res.mat) <- rownames(res.all)
  non.sig.count <- nrow(data) - sig.count
  rm(res.all)

  gc()
  RegisterData(dataSet)
  # record the sig gene vec
  return(c(1, sig.count, non.sig.count))
}

# perfor differential analysis for array/RNA seq data
# for two groups only (used for meta-analysis)
PerformLimma <- function(data, group) {
  library(limma)
  data <- data
  design <- model.matrix(~ -1 + group)
  fit <- lmFit(data, design)

  grps.cmp <- paste("group", levels(group)[2], " - ", "group", levels(group)[1], sep = "")
  myargs <- list(grps.cmp, levels = design)
  contrast.matrix <- do.call(makeContrasts, myargs)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  gc()
  return(list(fit.obj = fit))
}

# get result table from eBayes fit object
GetLimmaResTable <- function(fit.obj) {
  resTable <- topTable(fit.obj, number = Inf, adjust.method = "BH")
  if (!is.null(resTable$ID)) { # for older version
    rownames(resTable) <- resTable$ID
    resTable$ID <- NULL
  }
  return(resTable)
}

# given a gene id, plot its expression profile as box plot
PlotSelectedGene <- function(gene.id, type) {
  imgName <- paste("Gene_", gene.id, ".png", sep = "")
  library(lattice)
  if (anal.type == "onedata") {
    ids <- rownames(dataSet$resTable)
    inx <- which(ids == gene.id)
    symb <- dataSet$sig.genes.symbols[inx]
    if (type == "volcano") {
      symb <- ""
    }
    if (dataSet$comp.type == "custom") {
      Cairo(file = imgName, width = 280, height = 320, type = "png", bg = "white")
      grp.nms <- dataSet$grp.nms
      inx <- dataSet$cls %in% grp.nms
      cls <- dataSet$cls[inx]
      dat <- dataSet$data.norm[, inx]
      myplot <- bwplot(dat[gene.id, ] ~ as.character(cls),
        fill = "#0000ff22", scales = list(x = list(rot = 30)),
        xlab = "Class", ylab = "Expression Pattern", main = symb
      )
    } else if (length(dataSet$sec.cls) > 1) {
      out.fac <- as.character(dataSet$sec.cls)
      in.fac <- as.character(dataSet$fst.cls)
      colnames(dataSet$meta.info[, 1])

      Cairo(file = imgName, dpi = 72, width = 320, height = 320, type = "png", bg = "white")
      # ylim.ext <- GetExtendRange(dataSet$data.norm[gene.id, ], 12);
      layout <- c(2, 1)
      myplot <- bwplot(dataSet$data.norm[gene.id, ] ~ in.fac | out.fac,
        xlab = "Factors", ylab = "Expression Pattern", main = symb, scales = list(x = list(rot = 30)),
        fill = "#0000ff22", layout = layout
      )
    } else {
      Cairo(file = imgName, width = 280, height = 320, type = "png", bg = "white")

      myplot <- bwplot(dataSet$data.norm[gene.id, ] ~ as.character(dataSet$cls),
        fill = "#0000ff22", scales = list(x = list(rot = 30)),
        xlab = "Class", ylab = "Expression Pattern", main = symb
      )
    }
  } else { # metadata

    inmex.meta <- readRDS("inmex_meta.rds")
    if (inmex.meta$id.type == "entrez") {
      symb <- inmex.meta$gene.symbls[gene.id]
    } else {
      symb <- gene.id
    }
    num <- sum(mdata.all == 1)
    # calculate width based on the dateset number
    if (num == 1) {
      Cairo(file = imgName, width = 280, height = 320, type = "png", bg = "white")
      myplot <- bwplot(inmex.meta$plot.data[gene.id, ] ~ as.character(inmex.meta$cls.lbl),
        fill = "#0000ff22",
        xlab = "Class", ylab = "Expression Pattern", main = symb, scales = list(x = list(rot = 30))
      )
    } else {
      # calculate layout
      if (num < 6) {
        layout <- c(num, 1)
        height <- 320
        width <- 160 * num
      } else {
        rn <- round(num / 2)
        layout <- c(rn, 2)
        height <- 500
        width <- 160 * rn
      }

      Cairo(file = imgName, width = width, height = height, type = "png", bg = "white")
      data.lbl <- as.character(inmex.meta$data.lbl)
      data.lbl <- substr(data.lbl, 0, nchar(data.lbl) - 4)

      # get counts in each data, same order as a levels
      counts <- table(data.lbl)
      # back to factor
      data.lbl <- factor(data.lbl)

      # get new lbls to cut potential long names, and add sample numbers
      nlbls <- data.lbl
      levels(nlbls) <- abbreviate(levels(nlbls), 9)
      nlbls <- paste(levels(nlbls), "( n=", as.vector(counts), ")")
      # update labels
      data.lbl <- factor(data.lbl, labels = nlbls)
      # some time the transformed plot.data can switch class label, use the original data, need to be similar scale
      myplot <- bwplot(inmex.meta$plot.data[gene.id, ] ~ as.character(inmex.meta$cls.lbl) | data.lbl,
        xlab = "Datasets", ylab = "Expression Pattern", main = symb, scales = list(x = list(rot = 30)),
        fill = "#0000ff22", layout = layout
      )
    }
  }
  print(myplot)
  dev.off()
}

PlotSelectedGeneLoading <- function(gene.id) {
  if (anal.type == "metadata") {
    PlotSelectedGeneMeta(gene.id)
  } else {
    PlotSelectedGene(gene.id, "notVolcano")
  }
}

PlotSelectedGeneMeta <- function(gene.id) {

  # first get gene symbol
  inmex.meta <- readRDS("inmex_meta.rds")
  if (inmex.meta$id.type == "entrez") {
    symb <- inmex.meta$gene.symbls[gene.id]
  } else {
    symb <- gene.id
  }

  imgName <- paste("Gene_", gene.id, ".png", sep = "")
  library(lattice)

  num <- sum(mdata.all == 1)
  # calculate width based on the dateset number
  if (num == 1) {
    Cairo(file = imgName, width = 280, height = 320, type = "png", bg = "white")
    myplot <- bwplot(inmex.meta$plot.data[gene.id, ] ~ as.character(inmex.meta$cls.lbl),
      fill = "#0000ff22",
      xlab = "Class", ylab = "Expression Pattern", main = symb, scales = list(x = list(rot = 30))
    )
  } else {
    # this is a single long list
    layout <- c(1, num)
    height <- 200 * num
    width <- 280

    Cairo(file = imgName, width = width, height = height, type = "png", bg = "white")
    data.lbl <- as.character(inmex.meta$data.lbl)
    data.lbl <- substr(data.lbl, 0, nchar(data.lbl) - 4)

    # get counts in each data, same order as a levels
    counts <- table(data.lbl)
    # back to factor
    data.lbl <- factor(data.lbl)

    # get new lbls to cut potential long names, and add sample numbers
    nlbls <- data.lbl
    levels(nlbls) <- abbreviate(levels(nlbls), 9)
    nlbls <- paste(levels(nlbls), "( n=", as.vector(counts), ")")
    # update labels
    data.lbl <- factor(data.lbl, labels = nlbls)
    # some time the transformed plot.data can switch class label, use the original data, need to be similar scale
    myplot <- bwplot(inmex.meta$plot.data[gene.id, ] ~ as.character(inmex.meta$cls.lbl) | data.lbl,
      xlab = "Datasets", ylab = "Expression Pattern", main = symb, scales = list(x = list(rot = 30)),
      fill = "#0000ff22", layout = layout
    )
  }

  print(myplot)
  dev.off()
}

PlotCmpdView <- function(cmpdNm, format = "png", dpi = 72, width = NA) {
  if (anal.type == "onedata") {
    datanorm <- dataSet$data.norm
  } else {
    datanorm <- dataSet$data
  }
  clslbl <- dataSet$meta.info[, 1]
  imgName <- gsub("\\/", "_", cmpdNm)
  imgName <- paste(imgName, "_dpi", dpi, ".", format, sep = "")
  # indx<-which(rownames(boxplot_id)==cmpdNm);
  # gene.id <- boxplot_id[indx,1];
  gene.symb <<- doEntrez2SymbolMapping(cmpdNm)
  Cairo(file = imgName, dpi = dpi, width = 230, height = 230, type = format, bg = "transparent")
  par(mar = c(4, 3, 1, 2), oma = c(0, 0, 1, 0))
  boxplot(datanorm[which(rownames(datanorm) == as.character(cmpdNm)), ] ~ clslbl, las = 2, col = unique(GetColorSchema(clslbl)))
  title(main = gene.symb, out = T)
  dev.off()
  return(imgName)
}

# retrun the json obj
SaveHeatmapJSON <- function(fileName) {
  if (anal.type == "metadata") {
    json.res <- PrepareMetaHeatmapJSON()
  } else {
    json.res <- PrepareExpressHeatmapJSON()
  }

  library(RJSONIO)
  json.mat <- toJSON(json.res, .na = "null")
  sink(fileName)
  cat(json.mat)
  sink()
  current.msg <<- "Data is now ready for heatmap visualization!"
  return(1)
}


# prepare data for heatmap plotting include
# 1. meta info (top bars)
# 2. expression matrix
# 3. function annotation (left bars)
# all matrix will be combined,
# 1 and 2 separated by a row of 'null'
# 3 and 1+2 separated by a column of 'null'
PrepareExpressHeatmapJSON <- function() {
  sig.ids <- rownames(dataSet$sig.mat)
  stat.pvals <- dataSet$sig.mat$adj.P.Val

  # scale each gene
  data.stat <- readRDS("data.stat")
  hit.inz <- sig.ids %in% rownames(data.stat)
  sig.ids <- sig.ids[hit.inz]
  dat <- t(scale(t(data.stat[sig.ids, , drop = F])))

  # now pearson and euclidean will be the same after scaleing
  dat.dist <- dist(dat)

  orig.smpl.nms <- colnames(dat)
  orig.gene.nms <- rownames(dat)

  # do clustering and save cluster info
  # convert order to rank (score that can used to sort)
  if (nrow(dat) > 1) {
    dat.dist <- dist(dat)
    gene.ward.ord <- hclust(dat.dist, "ward.D")$order
    gene.ward.rk <- match(orig.gene.nms, orig.gene.nms[gene.ward.ord])
    gene.ave.ord <- hclust(dat.dist, "ave")$order
    gene.ave.rk <- match(orig.gene.nms, orig.gene.nms[gene.ave.ord])
    gene.single.ord <- hclust(dat.dist, "single")$order
    gene.single.rk <- match(orig.gene.nms, orig.gene.nms[gene.single.ord])
    gene.complete.ord <- hclust(dat.dist, "complete")$order
    gene.complete.rk <- match(orig.gene.nms, orig.gene.nms[gene.complete.ord])

    dat.dist <- dist(t(dat))
    smpl.ward.ord <- hclust(dat.dist, "ward.D")$order
    smpl.ward.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ward.ord])
    smpl.ave.ord <- hclust(dat.dist, "ave")$order
    smpl.ave.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ave.ord])
    smpl.single.ord <- hclust(dat.dist, "single")$order
    smpl.single.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.single.ord])
    smpl.complete.ord <- hclust(dat.dist, "complete")$order
    smpl.complete.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.complete.ord])
  } else {
    # force not to be single element vector which will be scaler
    stat.pvals <- matrix(stat.pvals)
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1)
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat)
  }

  gene.cluster <- list(
    pval = stat.pvals,
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk
  )

  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  )

  # prepare meta info
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  meta <- data.frame(dataSet$meta.stat)
  grps <- colnames(meta)
  nmeta <- meta.vec <- NULL
  0
  for (i in 1:ncol(meta)) {
    cls <- meta[, i]
    grp.nm <- grps[i]
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data
    ncls <- paste(grp.nm, as.numeric(cls)) # note, here to retain ordered factor
    nmeta <- c(nmeta, ncls)
  }

  # convert back to numeric
  nmeta <- as.numeric(as.factor(nmeta)) + 99
  unik.inx <- !duplicated(nmeta)

  # get corresponding names
  meta_anot <- meta.vec[unik.inx]
  names(meta_anot) <- nmeta[unik.inx] # name annotatation by their numbers

  nmeta <- matrix(nmeta, ncol = ncol(meta), byrow = F)
  colnames(nmeta) <- grps

  # for each gene/row, first normalize and then tranform real values to 30 breaks
  res <- t(apply(dat, 1, function(x) {
    as.numeric(cut(x, breaks = 30))
  }))

  # note, use {} will lose order; use [[],[]] to retain the order


  if (dataSet$annotated) {
    anot.id <- rownames(res)
    anot.res <- doEntrezIDAnot(anot.id)
    # single element vector will be converted to scalar, not array, need to prevent that
    gene.id <- anot.res$symbol
    if (length(gene.id) == 1) {
      gene.id <- matrix(gene.id)
    }
    gene.entrez <- anot.res$gene_id
    if (length(gene.entrez) == 1) {
      gene.entrez <- matrix(gene.entrez)
    }
    gene.name <- anot.res$name
    if (length(gene.name) == 1) {
      gene.name <- matrix(gene.name)
    }

    json.res <- list(
      data.type = dataSet$type,
      gene.id = anot.res$symbol,
      gene.entrez = gene.entrez,
      gene.name = anot.res$name,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    )
  } else if (file.exists("annotation.rds")) {
    # special gene.id and new gene.symbol
    anot.id <- rownames(res)
    anot.res <- doEntrezIDAnot(anot.id)
    gene.id <- rownames(anot.res)
    if (length(gene.id) == 1) {
      gene.id <- matrix(gene.id)
    }
    gene.entrez <- anot.res$gene_id
    if (length(gene.entrez) == 1) {
      gene.entrez <- matrix(gene.entrez)
    }
    gene.name <- paste(anot.res$symbol, anot.res$name, sep = " | ")
    if (length(gene.name) == 1) {
      gene.name <- matrix(gene.name)
    }

    json.res <- list(
      data.type = dataSet$type,
      gene.id = gene.id,
      gene.entrez = gene.entrez,
      gene.name = gene.name,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    )
  } else {
    gene.id <- orig.gene.nms
    if (length(gene.id) == 1) {
      gene.id <- matrix(gene.id)
    }
    json.res <- list(
      data.type = dataSet$type,
      gene.id = gene.id,
      gene.entrez = gene.id,
      gene.name = gene.id,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    )
  }
  return(json.res)
}

SaveClusterJSON <- function(fileNm, clustOpt, opt) {
  if (anal.type == "onedata") {
    SaveExpressClusterJSON(fileNm, clustOpt, opt)
  } else {
    initmetaloading <<- TRUE
    SaveMetaClusterJSON(fileNm, clustOpt, opt)
  }
}

SaveClusterJSONLoading <- function(fileNm, clustOpt, nb) {
  if (anal.type == "onedata") {
    SaveExpressClusterLoadingJSON(fileNm, clustOpt, nb)
  } else {
    SaveMetaClusterLoadingJSON(fileNm, clustOpt, nb)
  }
}

SaveExpressClusterLoadingJSON <- function(fileName, clustOpt, nb) {
  dat <- dataSet$data.norm
  pca3d <- list()
  dat <- na.omit(dat)
  nb <- as.numeric(nb)
  if (clustOpt == "pca") {
    pca <- prcomp(t(dat), center = T, scale = T)
    imp.pca <- summary(pca)$importance
    pca3d$score$axis <- paste("PC", 1:3, " (", 100 * round(imp.pca[2, ][1:3], 3), "%)", sep = "")
    coords <- data.frame(t(signif(pca$rotation[, 1:3], 5)))

    colnames(coords) <- NULL
    pca3d$score$xyz <- coords
    pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation))
    pca3d$score$entrez <- rownames(pca$rotation)
    weights <- imp.pca[2, ][1:3]
    mypos <- t(coords)
    meanpos <- apply(abs(mypos), 1, function(x) {
      weighted.mean(x, weights)
    })
    df <- data.frame(pos = meanpos, inx = seq.int(1, length(meanpos)))
    df <- df[order(-df$pos), ]
    if (nrow(df) > 2000) {
      inx <- df$inx[c(1:nb)]
      mypos <- mypos[inx, ]
      pca3d$score$xyz <- coords[inx]
      pca3d$score$name <- pca3d$score$name[inx]
      pca3d$score$entrez <- pca3d$score$entrez[inx]
    }
  }

  pca3d$cls <- dataSet$meta.info
  colnames(mypos) <- paste("Dim", 1:3, sep = "")
  # see if there is secondary
  loadEntrez <<- pca3d$score$entrez
  rownames(mypos) <- pca3d$score$name

  write.csv(mypos, file = "networkanalyst_3d_load_pos.csv")
  library(RJSONIO)
  json.mat <- toJSON(pca3d, .na = "null")
  sink(fileName)
  cat(json.mat)
  sink()
  current.msg <<- "Annotated data is now ready for PCA 3D visualization!"
  return(1)
}

# single expression data
SaveExpressClusterJSON <- function(fileName, clustOpt, opt) {
  dat <- dataSet$data.norm
  pca3d <- list()
  dat <- na.omit(dat)

  if (clustOpt == "pca") {
    if (opt == "all") {
      pca <- prcomp(t(dat), center = T, scale = T)
    } else {
      dat <- dat[which(rownames(dat) %in% loadEntrez), ]
      pca <- prcomp(t(dat), center = T, scale = T)
    }
    imp.pca <- summary(pca)$importance
    pca3d$score$axis <- paste("PC", 1:3, " (", 100 * round(imp.pca[2, ][1:3], 3), "%)", sep = "")
    coords <- data.frame(t(signif(pca$x[, 1:3], 5)))
  } else { # tsne
    require("Rtsne")
    dat <- as.matrix(t(dat))
    max.perx <- floor((nrow(dat) - 1) / 3)
    if (max.perx > 30) {
      max.perx <- 30
    }
    res <- Rtsne(dat, dims = 3, perplexity = max.perx)
    pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep = "")
    coords <- data.frame(t(signif(res$Y, 5)))
  }

  colnames(coords) <- NULL
  pca3d$score$xyz <- coords
  pca3d$score$name <- colnames(dataSet$data.norm)

  facA <- as.character(dataSet$fst.cls)
  if (all.numeric(facA)) {
    facA <- paste("Group", facA)
  }
  pca3d$score$facA <- facA

  mypos <- t(coords)
  colnames(mypos) <- paste("Dim", 1:3, sep = "")
  # see if there is secondary
  if (length(dataSet$sec.cls) > 1) {
    facB <- as.character(dataSet$sec.cls)
    if (all.numeric(facB)) {
      facB <- paste("Group", facB)
    }
    pca3d$score$facB <- facB

    # set shape based on the first group
    pca3d$score$shapes <- c("sphere", "triangle")

    # now set color based on 2nd group
    cols <- unique(GetColorSchema(dataSet$sec.cls))
    rgbcols <- col2rgb(cols)
    cols <- apply(rgbcols, 2, function(x) {
      paste("rgb(", paste(x, collapse = ","), ")", sep = "")
    })
    pca3d$score$colors <- cols

    mypos <- data.frame(factorA = facA, factorB = facB, mypos)
  } else {
    # now set color based on first group
    cols <- unique(GetColorSchema(dataSet$fst.cls))
    rgbcols <- col2rgb(cols)
    cols <- apply(rgbcols, 2, function(x) {
      paste("rgba(", paste(x, collapse = ","), ",1)", sep = "")
    })
    pca3d$score$colors <- cols
    mypos <- data.frame(factorA = facA, mypos)
  }
  pca3d$cls <- dataSet$meta.info
  rownames(mypos) <- colnames(dataSet$data.norm)

  write.csv(mypos, file = "networkanalyst_3d_pos.csv")
  library(RJSONIO)
  json.mat <- toJSON(pca3d, .na = "null")
  sink(fileName)
  cat(json.mat)
  sink()
  current.msg <<- "Annotated data is now ready for PCA 3D visualization!"
  return(1)
}

SetVolcanoHigh <- function(ids) {
  idsu <<- ids
  gene.vec <- unlist(strsplit(ids, "; "))
  volcanoHlVec <<- gene.vec
  if (length(volcanoHlVec) > 0) {
    return(1)
  } else {
    return(0)
  }
}

# prepare seeds from metaanalysis result
.prepareExpressSeeds <- function() {
  gene.list <- list()
  if (anal.type == "metadata") {
    if (selectedNetDataset == "meta_dat") {
      gene.mat <- meta.mat
      if (inmex.method != "votecount") {
        gene.list$metadata <- list(gene = rownames(gene.mat), logFC = unname(gene.mat[, 1]), adjP = unname(gene.mat[, 2]))
      } else {
        gene.list$metadata <- list(gene = rownames(gene.mat), logFC = unname(gene.mat[, 1]), adjP = unname(gene.mat[, 1]))
      }
    } else {
      # dataSet <- readRDS(selectedNetDataset);
      fit.obj.nm <- paste(selectedNetDataset, "fit.obj", sep = ".")
      fit2i <- readRDS(fit.obj.nm)
      gene.mat <- GetLimmaResTable(fit2i)
      gene.mat <- gene.mat[gene.mat$"adj.P.Val" < 0.05, ]
      gene.list$metadata <- list(gene = rownames(gene.mat), logFC = unname(gene.mat[, 1]), adjP = unname(gene.mat$"adj.P.Val"))
    }
  } else {
    gene.mat <- dataSet$sig.mat
    gene.list$metadata <- list(gene = rownames(gene.mat), logFC = unname(gene.mat$logFC), adjP = unname(gene.mat$"adj.P.Val"))
  }
  write.table(gene.mat, file = "sig_genes.txt")
  gene.vec <- rownames(gene.mat)
  GeneAnotDB <- doProteinIDMapping(rownames(gene.mat), "entrez")
  protein.vec <- GeneAnotDB[, 2]
  gene.mat <- data.matrix(gene.mat)
  rownames(gene.mat) <- protein.vec
  na.inx <- is.na(protein.vec)
  prot.mat <- gene.mat[!na.inx, , drop = F]
  # write.table(cbind(Uniprot=rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F);
  write.table(cbind(Emblprotein = rownames(prot.mat), Expression = prot.mat[, 1]), file = "seed_proteins.txt", row.names = F, quote = F)
  protein.vec <- prot.mat[, 1]
  if (length(protein.vec) == 1) {
    protein.vec <- as.matrix(protein.vec)
  }
  protein.list <- list()
  protein.list$metadata <- signif(protein.vec, 5)
  seed.expr <- prot.mat[, 1]
  seed.df <- as.matrix(seed.expr)
  rownames(seed.df) <- names(seed.expr)
  seed.expr <- RemoveDuplicates(seed.df, "max", quiet = F)
  seed.expr <<- seed.expr[, 1]
  seed.genes <<- unique(gene.vec)
  protein.vec <- unique(names(protein.vec))

  list(
    gene.list = gene.list,
    protein.list = protein.list,
    protein.vec = protein.vec
  )
}


# read tab delimited file
# can have many classes, stored in meta.info (starts with #)
# return a list (data.name, data.frame, meta.data)
.readTabData <- function(dataName) {
  if (length(grep("\\.zip$", dataName, perl = TRUE)) > 0) {
    dataName <- unzip(dataName)
    if (length(dataName) > 1) {
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep("MACOSX", dataName, perl = TRUE)
      if (length(osInx) > 0) {
        dataName <- dataName[-osInx]
      }
      dsInx <- grep("DS_Store", dataName, perl = TRUE)
      if (length(dsInx) > 0) {
        dataName <- dataName[-dsInx]
      }
      dat.inx <- grep(".[Tt][Xx][Tt]$", dataName)
      if (length(dat.inx) != 1) {
        current.msg <<- "More than one text files (.txt) found in the zip file."
        return(0)
      }
    }
  }

  NULL
  # using the powerful fread function, 10 times faster, note: default return data.table, turn off
  dat1 <- .readDataTable(dataName)

  # look for #CLASS, could have more than 1 class labels, store in a list
  meta.info <- list()
  cls.inx <- grep("^#CLASS", dat1[, 1])
  if (length(cls.inx) > 0) {
    for (i in 1:length(cls.inx)) {
      inx <- cls.inx[i]
      cls.nm <- substring(dat1[inx, 1], 2) # discard the first char #
      if (nchar(cls.nm) > 6) {
        cls.nm <- substring(cls.nm, 7) # remove class
      }
      cls.lbls <- dat1[inx, -1]
      # test NA
      na.inx <- is.na(cls.lbls)
      cls.lbls[na.inx] <- "NA"
      cls.lbls <- ClearFactorStrings(cls.nm, cls.lbls)

      meta.info[[cls.nm]] <- cls.lbls
    }
  } else {
    current.msg <<- "No metadata labels #CLASS found in your data!"
    return("F")
  }

  meta.info <- data.frame(meta.info)

  # now remove all comments in dat1
  # assign rownames after covert to matrix as data.frame does not allow duplicate names
  comments.inx <- grep("^#", dat1[, 1])
  dat1.nms <- dat1[-comments.inx, 1]
  dat1 <- dat1[-comments.inx, -1]
  dat1 <- data.matrix(dat1)
  rownames(dat1) <- dat1.nms

  list(
    name = basename(dataName),
    data = dat1,
    meta.info = meta.info
  )
}


# note, try to use the fread, however, it has issues with
# some windows 10 files "Line ending is \r\r\n. .... appears to add the extra \r in text mode on Windows"
# in such as, use the slower read.table method
.readDataTable <- function(fileName) {
  if (length(grep("\\.zip$", fileName, perl = TRUE)) > 0) {
    fileName <- unzip(fileName)
    if (length(fileName) > 1) {
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep("MACOSX", fileName, perl = TRUE)
      if (length(osInx) > 0) {
        fileName <- fileName[-osInx]
      }
      dsInx <- grep("DS_Store", fileName, perl = TRUE)
      if (length(dsInx) > 0) {
        fileName <- fileName[-dsInx]
      }
      dat.inx <- grep(".[Tt][Xx][Tt]$", fileName)
      if (length(dat.inx) != 1) {
        current.msg <<- "More than one text files (.txt) found in the zip file."
        return(0)
      }
    }
  }
  dat <- try(data.table::fread(fileName, header = TRUE, check.names = FALSE, data.table = FALSE))
  if (class(dat) == "try-error") {
    # try to use "tr" to remove double return characters
    trFileName <- paste("tr -d \'\\r\' <", fileName)
    dat <- try(data.table::fread(trFileName, header = TRUE, check.names = FALSE, data.table = FALSE))
    if (class(dat) == "try-error") {
      print("Using slower file reader ...")
      formatStr <- substr(fileName, nchar(fileName) - 2, nchar(fileName))
      if (formatStr == "txt") {
        dat <- try(read.table(fileName, header = TRUE, comment.char = "", check.names = F, as.is = T))
      } else { # note, read.csv is more than read.table with sep=","
        dat <- try(read.csv(fileName, header = TRUE, comment.char = "", check.names = F, as.is = T))
      }
    }
  }
  return(dat)
}

meanSdPlot <- function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
                       ylab = "sd", pch, plot = TRUE, bins = 50, ...) {
  stopifnot(is.logical(ranks), length(ranks) == 1, !is.na(ranks))

  n <- nrow(x)
  if (n == 0L) {
    warning("In 'meanSdPlot': input matrix 'x' has 0 rows. There is nothing to be done.")
    return()
  }
  if (!missing(pch)) {
    warning("In 'meanSdPlot': 'pch' is ignored.")
  }

  px <- rowMeans(x, na.rm = TRUE)
  py <- sqrt(rowV(x, mean = px, na.rm = TRUE))
  rpx <- rank(px, na.last = FALSE, ties.method = "random")

  ## run median with centers at dm, 2*dm, 3*dm,... and width 2*dm
  dm <- 0.025
  midpoints <- seq(dm, 1 - dm, by = dm)
  within <- function(x, x1, x2) {
    (x >= x1) & (x <= x2)
  }
  mediwind <- function(mp) median(py[within(rpx / n, mp - 2 * dm, mp + 2 * dm)], na.rm = TRUE)
  rq.sds <- sapply(midpoints, mediwind)

  res <- if (ranks) {
    list(rank = midpoints * n, sd = rq.sds, px = rpx, py = py)
  } else {
    list(quantile = quantile(px, probs = midpoints, na.rm = TRUE), sd = rq.sds, px = px, py = py)
  }

  fmt <- function() function(x) format(round(x, 0), nsmall = 0L, scientific = FALSE)

  res$gg <- ggplot(
    data.frame(px = res$px, py = res$py),
    aes_string(x = "px", y = "py")
  ) + xlab(xlab) + ylab(ylab) +
    geom_hex(bins = bins, ...) +
    scale_fill_gradient(name = "count", trans = "log", labels = fmt()) +
    geom_line(aes_string(x = "x", y = "y"),
      data = data.frame(x = res[[1]], y = res$sd), color = "red"
    )

  if (plot) print(res$gg)

  return(invisible(res))
}

rowV <- function(x, mean, ...) {
  sqr <- function(x) x * x
  n <- rowSums(!is.na(x))
  n[n < 1] <- NA
  if (missing(mean)) {
    mean <- rowMeans(x, ...)
  }
  return(rowSums(sqr(x - mean), ...) / (n - 1))
}

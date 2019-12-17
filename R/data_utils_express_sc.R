
##################################################
## R script for NetworkAnalyst
## Description: functions for single-cell RNA-Seq expression data
##
## Author: Matthew J. Oldach, matthew-oldach@mcgill.ca
###################################################

# read tab delimeted file
# convert to Seurat object which
# contains a matrix of counts and metadata
ReadTabExpressData_sc <- function(fileName) {
  if (missing(fileName)) {
    msg <- paste0("The required count table was not provided. Please select a valid file.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } # was the required argument supplied?

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
    msg <- c(msg, paste("removed 5% features with near-constant values"))
  }

  minVal <- min(int.mat, na.rm = T)
  na.inx <- is.na(int.mat)
  if (sum(na.inx) > 0) {
    int.mat[na.inx] <- minVal / 2
    msg <- c(msg, "the remaining", sum(na.inx), "missing variables were replaced with data min")
  }
  current.msg <<- paste(msg, collapse = "; ")
  data.proc <- RemoveDuplicates(int.mat, "mean", quiet = T) # WHICH PACKAGE DOES THIS FUNCTION COME FROM?!?

  # save processed data for download user option
  write.csv(data.proc, file = "data_processed.csv")
  saveRDS(data.proc, "data.proc.rds")
  dataSet <<- dataSet
  return(1)
}

PerformDataAnnot_sc <- function(org, dataType, idType, lvlOpt) {
  if (missing(org)) {
    msg <- paste0("Please provide a valid species for the 'org' parameter.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(dataType)) {
    msg <- paste0("Please provide a valid 'dataType' option.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(idType)) {
    msg <- paste0("Please provide a valid 'idType' option.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(lvlOpt)) {
    msg <- paste0("Please provide a valid 'lvlOpt' option.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }

  data.org <<- org
  SetInitLib(org)

  dataSet$type <- dataType
  dataSet$id.orig <- dataSet$id.current <- idType
  dataSet$annotated <- F
  # should not contain duplicates, however sanity check
  data.proc <- as.data.frame(readRDS("data.proc.rds"))
  dataSet$data.anot <- data.proc

  # get counts
  totalCount <- sum(colSums(data.proc))
  avgCount <- sum(colSums(data.proc)) / ncol(data.proc)
  minCount <- min(colSums(data.proc))
  maxCount <- max(colSums(data.proc))

  if (org != "NA" & idType != "NA") {
    feature.vec <- rownames(data.proc)
    anot.id <- doAnnotation(feature.vec, idType)

    # dataSet$annotation <- anot.id;
    saveRDS(anot.id, "annotation.rds")

    hit.inx <- !is.na(anot.id)
    matched.len <- sum(hit.inx)
    perct <- round(matched.len / length(feature.vec), 3) * 100
    thresh <- 0.1 # previous value of 0.25 is causing challenges
    # if no entrez match, keep ensemblgene_id
    na.inx <- is.na(anot.id)
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
        data.proc[hit.inx, ] -> matched_df

        unmatched.entrez <- anot.id[na.inx]
        data.proc[na.inx, ] -> unmatched_df

        data.anot <- rbind(unmatched_df, matched_df)

        ## need a step to remove duplicate entrez ID's
        data.anot$entrez <- c(matched.entrez, unmatched.entrez)

        if (lvlOpt == "mean") {
          data.anot <- aggregate(data.anot[, -which(names(data.anot) %in% "entrez")], by = list(entrez = data.anot$entrez), FUN = mean)
        } else {
          data.anot <- aggregate(data.anot[, -which(names(data.anot) %in% "entrez")], by = list(entrez = data.anot$entrez), FUN = sum)
        }
        rownames(data.anot) <- data.anot$entrez
        data.anot <- data.anot[, -which(names(data.anot) %in% "entrez")]

        current.msg <<- paste(current.msg, "Data is now transformed to gene-level (Entrez) expression.")
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
    current.msg <<- paste("No annotation was performed. Make sure organism and gene ID are specified correctly!")
  }
  # need to save the ids (mixed gene annotation and original id)
  # in case, users needs to keep unannotated features
  # this need to be updated together with data from now on
  dataSet <<- dataSet

  saveRDS(data.anot, file = "annotation.rds")


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

MakeSeuratObject <- function(min.cells = min.cells, min.features = min.features, filterUnmapped = true) {
  min.cells <- as.numeric(min.cells)
  min.features <- as.numeric(min.features)

  if (missing(min.cells)) {
    msg <- paste0("Please provide a valid 'min.cells' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(min.features)) {
    msg <- paste0("Please provide a valid 'min.features' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  if (!is.numeric(min.cells)) {
    msg <- paste0("'min.cells' must be numeric.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (!is.numeric(min.features)) {
    msg <- paste0("'min.features' must be numeric.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (min.cells < 0) {
    msg <- paste0("'min.cells' must be positive.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (min.features < 0) {
    msg <- paste0("'min.features' must be positive.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }

  library(Seurat)
  if (filterUnmapped == "false") {
    dataSet <- readRDS("data.proc.rds")
  } else {
    dataSet <- readRDS("annotation.rds")
  }

  dataSet <- CreateSeuratObject(counts = dataSet, min.cells = min.cells, min.features = min.features) # add nCount_RNA & nFeature_RNA

  dataSet[["percent.mt"]] <- PercentageFeatureSet(dataSet, pattern = "^MT-") # add percentage mitochondria - remove this!

  saveRDS(dataSet, "data.seurat.rds")
  return(1)
}

# show counts vs gene per cell
PlotCellSummary <- function(imgNm, dpi, format) {
  if (missing(imgNm)) {
    msg <- paste0("Please provide a valid 'imgNm' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")


  counts <- data.table::fread("data_processed.csv")

  col_2_row <- function(.data, var = "rowname") {
    .data <- as.data.frame(.data)
    rownames(.data) <- .data[[var]]
    .data[[var]] <- NULL
    .data
  }

  counts <- col_2_row(counts, var = "V1")

  # cell summaries
  plot_cell <- data.frame(
    counts_per_cell = colSums(counts),
    genes_per_cell = colSums(counts > 0)
  )

  Cairo(file = imgNm, width = 600 * dpi / 72, height = 450 * dpi / 72, type = format, bg = "white", dpi = dpi, unit = "px")
  plot(plot_cell$counts_per_cell, plot_cell$genes_per_cell, log = "xy", col = "#F8766D", pch = 16)
  title("counts vs genes per cell")

  dev.off()
}

# visualize QC metrics, and use these to filter cells
PlotQC <- function(imgNm, dpi, format) {
  if (missing(imgNm)) {
    msg <- paste0("Please provide a valid 'imgNm' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")

  library("Seurat")
  dataSet <- readRDS("data.seurat.rds")

  Cairo(file = imgNm, width = 600 * dpi / 72, height = 450 * dpi / 72, type = format, bg = "white", dpi = dpi, unit = "px")
  g <- VlnPlot(dataSet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(g)
  dev.off()
}

# visualize feature-feature relationships
PlotScatter <- function(imgNm, dpi, format) {
  if (missing(imgNm)) {
    msg <- paste0("Please provide a valid 'imgNm' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")

  library("Seurat")

  dataSet <- readRDS("data.seurat.rds")

  Cairo(file = imgNm, width = 600 * dpi / 72, height = 450 * dpi / 72, type = format, bg = "white", dpi = dpi, unit = "px")
  plot1 <- FeatureScatter(dataSet, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(dataSet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  g <- CombinePlots(plots = list(plot1, plot2))
  print(g)
  dev.off()
}

PerformExpressNormalization_sc <- function(norm.opt = norm.opt, feature.min = 200, feature.max = 2500, mt.min = 5, scaleSC = TRUE, scale.factor = 10000, hvfs = 2000) {
  save.image("whats-wrong.rds")
  if (missing(norm.opt)) {
    msg <- paste0("Please provide a valid 'norm.opt' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(feature.min)) {
    msg <- paste0("Please provide a valid 'feature.min' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(feature.max)) {
    msg <- paste0("Please provide a valid 'feature.max' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(mt.min)) {
    msg <- paste0("Please provide a valid 'mt.min' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(scale.factor)) {
    msg <- paste0("Please provide a valid 'scale.factor' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(scale)) {
    msg <- paste0("Please provide a valid 'scale' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (missing(hvfs)) {
    msg <- paste0("Please provide a valid 'hvfs' value.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  if (!is.numeric(feature.min)) {
    msg <- paste0("'feature.min' must be numeric.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (!is.numeric(feature.max)) {
    msg <- paste0("'feature.max' must be numeric.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (!is.numeric(mt.min)) {
    msg <- paste0("'mt.min' must be numeric.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (!is.numeric(scale.factor)) {
    msg <- paste0("'scale.factor' must be numeric.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (!is.numeric(hvfs)) {
    msg <- paste0("'hvfs' must be numeric.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (!is.logical(scale)) {
    msg <- paste0("'hvfs' must be logical.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (!is.character(norm.opt)) {
    msg <- paste0("'norm.opt' must be a character")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (feature.min < 0) {
    msg <- paste0("'feature.min' must be zero or positive.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (feature.max < 0) {
    msg <- paste0("'feature.max' must be zero or positive.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (feature.min > feature.max) {
    msg <- paste0("'feature.min' must be less than 'feature.max'.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (mt.min > 100) {
    msg <- paste0("'mt.min' must be between 0 - 100%.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (mt.min < 0) {
    msg <- paste0("'mt.min' must be between 0 - 100%.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (scale.factor < 0) {
    msg <- paste0("'scale.factor' must be positive.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  } else if (hvfs < 0) {
    msg <- paste0("'hvfs' must be positive.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }

  library("Seurat")
  dataSet <- readRDS("data.seurat.rds")

  # dataSet <- Seurat:::subset.Seurat(dataSet, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # cannot hard-code: https://github.com/satijalab/seurat/issues/1053#issuecomment-454512002
  nFeature_RNA <- FetchData(object = dataSet, vars = "nFeature_RNA")
  dataSet <- dataSet[, which(x = nFeature_RNA > feature.min & nFeature_RNA < feature.max)]
  percent.mt <- FetchData(object = dataSet, vars = "percent.mt")
  dataSet <- dataSet[, which(x = percent.mt < mt.min)]

  # The default # npcs is 50, if there are fewer than 50 cells set to one less than
  number_of_cells <- ncol(dataSet)
  if (number_of_cells < 50) {
    max.npcs <- number_of_cells - 1
  } else {
    max.npcs <- 50
  }

  if (norm.opt != "LogNormalize" && norm.opt != "CLR" && norm.opt != "RC") {
    msg <- paste0("Please enter a valid normalization method.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dataSet <- Seurat::NormalizeData(dataSet, normalization.method = norm.opt, scale.factor = scale.factor) # two other normalization methods: CLR and RC

  # Seurat::FindVariableFeatures() has a number of parameters we may want to include in another function in the future...
  dataSet <- FindVariableFeatures(dataSet, selection.method = "vst", nfeatures = hvfs) # would users want a different selection.method?

  if (scale == TRUE) {
    all.genes <- rownames(dataSet)
    dataSet <- ScaleData(dataSet, features = all.genes)
    dataSet <- RunPCA(dataSet, features = VariableFeatures(object = dataSet), npcs = max.npcs)
  } else {
    dataSet <- RunPCA(dataSet, features = VariableFeatures(object = dataSet), npcs = max.npcs)
  }
  # save normalized data for download user option
  saveRDS(dataSet, "data.normalized.rds") # save Seurat object
  data.table::fwrite(as.data.frame(dataSet@assays[["RNA"]]@scale.data), file = "data_normalized.csv") # save count matrix
  return(1)
}

# elbow plot
PlotElbow <- function(imgNm, dpi, format) {
  if (missing(imgNm)) {
    msg <- paste0("Please provide a vaild 'imgNm'.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")

  library(Seurat)
  dataSet <- readRDS("data.normalized.rds")

  g <- ElbowPlot(dataSet)
  Cairo(file = imgNm, width = 600 * dpi / 72, height = 450 * dpi / 72, type = format, bg = "white", dpi = dpi, unit = "px")
  print(g)
  dev.off()
}

ChoosePCS <- function() {
  library(Seurat)
  dataSet <- readRDS("data.normalized.rds")
  if (scale == TRUE) {
    dataSet <- FindNeighbors(dataSet, reduction = "pca", dims = 1:20)
    dataSet <- FindClusters(dataSet, resolution = 0.5, algorithm = 1)
    # re-assign identities as cluster found by Seurat
    Idents(dataSet) <- dataSet@meta.data$seurat_clusters
  } else {
    # graph-based clustering approach
    dataSet <- FindNeighbors(dataSet, reduction = "pca", dims = 1:20)
    dataSet <- FindClusters(dataSet, resolution = 0.5, algorithm = 1)
    # re-assign identities as cluster found by Seurat
    Idents(dataSet) <- dataSet@meta.data$seurat_clusters
  }
  saveRDS(dataSet, "data.pcs.rds")
}

PlotTSNE <- function(imgNm, dpi, format) {
  if (missing(imgNm)) {
    msg <- paste0("Please provide a vaild 'imgNm'.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")

  library(Seurat)
  dataSet <- readRDS("data.pcs.rds")

  ## Rather than naively trying a loop that trys smaller values of perplexity (from 50 -> 5)
  ## lets write a function that calls the function we're interested in, that returns
  ## a numerical value indicating how close on is to the "best" result achievable
  ## (e.g. largest perplexity value that returns a computed result)
  f <- function(i) {
    zzz <- try(RunTSNE(object = dataSet, dims.use = 1:10, do.fast = TRUE, perplexity = i), TRUE)
    if (class(zzz) == "try-error") {
      "ok!"
    } else {
      stop("yuck")
    }
  }

  g <- function(i) {
    tryCatch(
      {
        f(i)
        i
      },
      error = function(e) {
        0
      }
    )
  }

  z <- sapply(70:5, g) # if you find this erroring-out improve code logic here
  z <- z[z > 0]
  h <- min(z) - 1
  dataSet <- RunTSNE(object = dataSet, dims.use = 1:10, do.fast = TRUE, perplexity = h)

  saveRDS(dataSet, "data.tSNE.rds")

  g <- DimPlot(object = dataSet, reduction = "tsne") # note that you can set do.label=T to help label individual clusters

  Cairo(file = imgNm, width = 600 * dpi / 72, height = 450 * dpi / 72, type = format, bg = "white", dpi = dpi, unit = "px")

  print(g)

  dev.off()
}

# umap of scaled data
PlotUmap <- function(imgNm, dpi, format) {
  if (missing(imgNm)) {
    msg <- paste0("Please provide a vaild 'imgNm'.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")

  library("Seurat")
  dataSet <- readRDS("data.tSNE.rds")

  # number of cells left
  number_of_cells <- ncol(dataSet)
  n_neighbors <- number_of_cells - 1

  dataSet <- RunUMAP(dataSet, reduction = "pca", dims = 1:20, n_neighbors = n_neighbors)
  saveRDS(dataSet, "data.tSNE.rds")

  g <- DimPlot(dataSet, reduction = "umap", split.by = "seurat_clusters")

  Cairo(file = imgNm, width = 600 * dpi / 72, height = 450 * dpi / 72, type = format, bg = "white", dpi = dpi, unit = "px")

  print(g)

  dev.off()
}

# calculate a subset of features that exhibit high cell-to-cell variation in the dataset
PlotHvfs <- function(imgNm, dpi, format) {
  if (missing(imgNm)) {
    msg <- paste0("Please provide a vaild 'imgNm'.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep = "")
  library("Seurat")

  dataSet <- readRDS("data.tSNE.rds")

  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(dataSet), 10)

  # plot variable features with labels
  plot1 <- VariableFeaturePlot(dataSet)
  g <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

  Cairo(file = imgNm, width = 600 * dpi / 72, height = 450 * dpi / 72, type = format, bg = "white", dpi = dpi, unit = "px")

  print(g)

  dev.off()
}

# note, setup the main class, keep the original order
SetMainClass_sc <- function() {
  dataSet <- readRDS("data.tSNE.rds")
  lbls <- factor(dataSet@meta.data[["seurat_clusters"]])
  lvls.orig <- unique(lbls)
  cls <- factor(lbls, levels = lvls.orig, ordered = T)
  cls # record main cls
  return(levels(cls))
}

# perform differential analysis
# default: compare one group against all other cells
# custom: only compare two groups (A-C)
PerformDEAnal_sc <- function(imgNm, anal.type = "default", par1 = NULL, par2 = NULL, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = FALSE) {
  library(dplyr)
  library(Seurat)
  # As a default, Seurat performs differential expression based on the non-parameteric Wilcoxon rank sum test
  dataSet <- readRDS("data.tSNE.rds")

  if (anal.type != "default" && anal.type != "cluster" && method != "RC") {
    msg <- paste0("Please provide a vaild 'anal.type'.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  if (test.use != "wilcox" && test.use != "bimod" && test.use != "t" && test.use != "negbinom" && test.use != "poisson" && test.use != "LR" && test.use != "MAST" && test.use != "DESeq2" && test.use != "roc") {
    msg <- paste0("Please provide a vaild DE test in 'test.use'.")
    print(msg)
    norm.msg <<- current.msg <<- msg
    return(0)
  }
  if (anal.type == "default") {
    # find topFeatures for every cluster compared to all remaining cells, report only the positive one
    topFeatures <- FindAllMarkers(dataSet, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold)
    colnames(topFeatures)[7] <- "entrezgene_id"
    # save the list of candidate marker genes for further examination
    data.table::fwrite(topFeatures, "candidate-marker-genes.txt")
    topFeatures %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_logFC)

    # add metadata
    filename <<- "SigGene_pairwise"
    resTable <<- topFeatures

    # make heatmap
    g <- DoHeatmap(dataSet, features = c(head(rownames(topFeatures), 10))) + NoLegend()
    png(file = paste0(imgNm, "_heatmap.png"), bg = "white", unit = "px")
    print(g)
    dev.off()
    # make violin plots
    g2 <- VlnPlot(dataSet, features = c(head(rownames(topFeatures), 3)))
    png(file = paste0(imgNm, "_violinplot.png"), bg = "white", unit = "px")
    print(g2)
    dev.off()
  } else if (anal.type == "cluster") {
    # find all topFeatures for cluster of interest
    topFeatures <- FindMarkers(dataSet, ident.1 = par1, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold)
    colnames(topFeatures)[7] <- "entrezgene_id"
    # save the list of candidate marker genes for further examination
    data.table::fwrite(topFeatures, "candidate-marker-genes.txt")
    topFeatures %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_logFC)

    # add metadata
    filename <<- paste0("SigGene_", par1, "_vs_Else")
    resTable <<- topFeatures

    g <- DoHeatmap(dataSet, features = c(head(rownames(topFeatures), 10))) + NoLegend()
    png(file = paste0(imgNm, "_heatmap.png"), bg = "white", unit = "px")
    print(g)
    dev.off()
    # make violin plots
    g2 <- VlnPlot(dataSet, features = c(head(rownames(topFeatures), 3)))
    png(file = paste0(imgNm, "_violinplot.png"), bg = "white", unit = "px")
    print(g2)
    dev.off()
  }
}

# update result based on new cutoff
GetSigGenes_sc <- function(res.nm, p.lvl, fc.lvl, update = T, inx) {
  resTable <- data.table::fread("candidate-marker-genes.txt")
  dataSet <- readRDS("data.tSNE.rds")
  total <- nrow(dataSet)
  filename <- filename

  if (nrow(resTable) == 0) {
    current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff.")
  }

  resTable <- resTable[ which(resTable$p_val_adj < p.lvl), ]

  de.Num <- nrow(resTable)

  # display at most 5000 genes for the server (two main reasons)
  # 1) should not have more 22% (human: 23000) DE of all genes (biological)
  # 2) IE canvas can display no more than 6800 pixels (computational)
  if (nrow(resTable) > 5000) {
    resTable <- resTable[1:5000, ]
    current.msg <<- paste(current.msg, " Due to computational constraints, only top 5000 genes will be used. ", collapse = "\n")
  }

  saveRDS(data, file = "data.stat")

  anot.id <- rownames(resTable)
  doEntrezIDAnot(anot.id)
  # save here!?

  gene <- rownames(resTable)
  logFC <- resTable$avg_logFC

  geneList <- paste(gene, logFC, collapse = "\n")
  resTable[ which(resTable$p_val_adj < p.lvl), ]
  up <- nrow(resTable[ which(resTable$avg_logFC > 2), ])
  down <- nrow(resTable[ which(resTable$avg_logFC < 2), ])

  dataSet <<- dataSet

  lst <- list(colnames(dataSet), dataSet@assays[["RNA"]]@scale.data, data.frame(CLASS = dataSet@meta.data[["seurat_clusters"]]), resTable, rownames(dataSet), org = data.org)

  library(RJSONIO)
  json.obj <- toJSON(lst)
  sink("NetworkAnalyst_matrix.json")
  cat(json.obj)
  return(c(filename, de.Num, geneList, total, up, down))
}

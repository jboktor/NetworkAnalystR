##################################################
## R script for NetworkAnalyst
## Description: prepare data for chord diagram
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

PrepareChordDataInit <- function() {
  # create a list store all possible combination (for a max of 4)
  # note, the whole region are divided into small areas (16 for 4; 7 for 3, 3 for 2)
  # for instance:
  # a: a unique (no b, no c)
  # ab: a and b, no c
  list()
  if (anal.type == "metadata") {
    sel.inx <- mdata.all == 1
    sel.nms <- names(mdata.all)[sel.inx]
  } else {
    sel.nms <- listNms
  }
  venn.genenb <- vector()
  sel.dats <- list()
  for (i in 1:length(sel.nms)) {
    nm <- sel.nms[i]
    dataSet <- readRDS(nm)
    if (anal.type == "metadata") {
      sel.dats[[nm]] <- rownames(dataSet$sig.mat)
    } else {
      sel.dats[[nm]] <- rownames(dataSet$prot.mat)
    }
    venn.genenb[i] <- length(sel.dats[[nm]])
  }
  if (anal.type == "metadata" & meta.selected) {
    sel.dats[["meta_dat"]] <- as.character(meta.stat$de)
    venn.genenb[length(venn.genenb) + 1] <- length(as.character(meta.stat$de))
  }

  chord.list <<- sel.dats
  chord.genenb <<- venn.genenb
  chord.list.up <<- sel.dats
  chord.genenb.up <<- venn.genenb
  return(PrepareChordData())
}

PrepareSelChordData <- function(selectedNms) {
  list()
  sel.nms <- unlist(strsplit(selectedNms, ";"))
  nm.vec <<- sel.nms
  SelectData()
  venn.genenb <- vector()
  sel.dats <- list()
  for (i in 1:length(sel.nms)) {
    nm <- sel.nms[i]
    if (nm != "meta_dat") {
      dataSet <- readRDS(nm)
      if (anal.type == "metadata") {
        sel.dats[[nm]] <- rownames(dataSet$sig.mat)
      } else {
        sel.dats[[nm]] <- rownames(dataSet$prot.mat)
      }
      venn.genenb[i] <- length(sel.dats[[nm]])
    } else {
      sel.dats[[nm]] <- as.character(meta.stat$de)
      venn.genenb[i] <- length(as.character(meta.stat$de))
    }
  }

  chord.list.up <<- sel.dats
  chord.genenb.up <<- venn.genenb
  return(PrepareChordData())
}

GetSelectedDataNamesChord <- function() {
  return(paste(names(chord.list), collapse = ";"))
}

GetSelectedDataGeneNumberChord <- function() {
  return(paste(chord.genenb, collapse = ";"))
}

GetSelectedDataNamesUpdatedChord <- function() {
  return(paste(names(chord.list.up), collapse = ";"))
}

GetSelectedDataGeneNumberUpdatedChord <- function() {
  return(paste(chord.genenb.up, collapse = ";"))
}

GetChordFileCount <- function() {
  return(chord_count - 1)
}

PrepareChordData <- function() {
  if (anal.type == "metadata") {
    res <- PrepareMetaChordData()
  } else {
    res <- PrepareListChordData()
  }

  if (is.null(res)) {
    return(0)
  }
  library(RJSONIO)
  chordData <- res$chordData
  fileNm <- paste0("networkanalyst_chorddata_", chord_count, ".json")
  sink(fileNm)
  cat(toJSON(chordData))
  sink()

  lookup <- res$lookup
  fileNm2 <- paste0("networkanalyst_chord_lookup_", chord_count, ".json")
  sink(fileNm2)
  cat(toJSON(lookup))
  chord_count <<- chord_count + 1
  sink()
  return(1)
}

PrepareChordDataFromList <- function(newDat, uniq.enIDs) {
  # now combine lists into a single matrix
  all.nms <- listNms
  nm.used <- NULL
  for (i in 1:length(all.nms)) {
    dataNm <- all.nms[i]
    if (mdata.all[[dataNm]] == 1) { # selected
      nm.used <- c(nm.used, dataNm)
    }
  }
  if (anal.type == "metadata" & meta.selected) {
    nm.used <- c(nm.used, "meta")
  }
  hit.mat <- fc.mat <- matrix(0, nrow = length(nm.used), ncol = length(uniq.enIDs))
  nm.mat <- matrix(NA, nrow = length(nm.used), ncol = length(uniq.enIDs))
  shared.ids <- uniq.enIDs
  shared.sbls <- doEntrez2SymbolMapping(uniq.enIDs)

  dataNms <- nm.used
  colnames(hit.mat) <- colnames(fc.mat) <- colnames(nm.mat) <- shared.sbls
  rownames(hit.mat) <- rownames(fc.mat) <- rownames(nm.mat) <- dataNms
  names(shared.sbls) <- shared.ids

  # push to parent env.
  shared.sbls <<- shared.sbls

  for (i in 1:length(dataNms)) {
    nm <- dataNms[i]
    exp.vec <- newDat[[nm]]
    hit.inx <- shared.ids %in% names(exp.vec)
    hit.mat[i, hit.inx] <- 1
    nm.mat[i, hit.inx] <- nm
    fc.mat[i, hit.inx] <- as.numeric(exp.vec)
  }

  ############## Links#######
  ## links & arcs data
  ###########################

  links <- arcs <- list()

  # arcs data
  de.num <- apply(hit.mat, 1, sum)
  de.prct <- as.numeric(de.num / sum(de.num))

  for (i in 1:length(dataNms)) {
    arcs[[i]] <- list(
      name = dataNms[i],
      label = dataNms[i],
      size = as.numeric(de.num[i]),
      val = de.prct[i]
    )
  }

  # now, get unique names by each inserting in the gene names
  lnk.genes <- colnames(nm.mat)
  nm.mat <- apply(
    nm.mat, 1,
    function(x) {
      gd.hits <- !is.na(x)
      x[gd.hits] <- paste(x[gd.hits], "*", lnk.genes[gd.hits], sep = "")
      x
    }
  )

  nm.mat <- unname(t(nm.mat))
  for (l in 1:nrow(nm.mat)) {
    x <- nm.mat[l, ]
    hit.inx <- !is.na(x)
    links$name <- c(links$name, x[hit.inx])

    y.mat <- nm.mat[, hit.inx, drop = F]
    new.mat <- apply(
      y.mat, 2,
      function(d) {
        myd <- d[-l]
        if (all(is.na(myd))) { # link to itself?
          d[l]
        } else {
          myd[!is.na(myd)]
        }
      }
    )
    links$imports <- c(links$imports, new.mat)
  }

  # need to reformat
  links.new <- list()
  for (l in 1:length(links$name)) {
    impt <- links$imports[[l]]
    if (length(impt) == 1) {
      impt <- as.matrix(impt)
    }
    links.new[[l]] <- list(
      name = links$name[l],
      imports = impt
    )
  }

  ###################
  # colors & labels
  library(RColorBrewer)
  if (length(dataNms) > 9) {
    colors <- brewer.pal(length(dataNms), "Set3")
  } else if (length(dataNms) > 2) {
    colors <- brewer.pal(length(dataNms), "Dark2")
  } else {
    colors <- c("#7570B3", "#D95F02")
  }
  labels <- rep("", length(dataNms))
  names(labels) <- names(colors) <- dataNms

  ###################
  # links done!, now get weights/pvalues/data for every pairs of chord
  dataNms <- rownames(fc.mat)
  geneNms <- colnames(fc.mat)
  entrezIDs <- names(shared.sbls)
  weights <- data <- list()

  # also need to setup lookup json
  lookup <- list()
  for (i in 1:length(dataNms)) {
    nlst <- vector(length = length(dataNms) - 1, mode = "list")
    nm <- dataNms[i]
    names(nlst) <- dataNms[-i]
    lookup[[nm]] <- nlst
  }

  item.count <- 0
  for (i in 1:ncol(nm.mat)) {
    geneNm <- geneNms[i]
    entrez <- entrezIDs[i]
    hit.inx <- !is.na(nm.mat[, i])
    fc.vals <- as.numeric(fc.mat[hit.inx, i])

    datNms <- dataNms[hit.inx]

    if (length(datNms) == 1) {
      item.count <- item.count + 1
      id <- paste(datNms, "*", datNms, "*", geneNm, sep = "")
      # weights[[id]]<-fc.vals;
      weights[[id]] <- 1
      lookup[[datNms]][[datNms]][[geneNm]] <- list(
        entrez = entrez,
        genesymbol = geneNm,
        weight = fc.vals
      )
      data[[item.count]] <- list(
        fromModule = datNms,
        toModule = datNms,
        arc = geneNm,
        weight = fc.vals,
        fromColor = as.character(colors[datNms]),
        toColor = as.character(colors[datNms])
      )
    } else {
      for (m in 1:(length(datNms) - 1)) {
        nm1 <- datNms[m]
        fc.vals[m]
        for (n in (m + 1):length(datNms)) {
          item.count <- item.count + 1
          nm2 <- datNms[n]
          fc.vals[n]
          id1 <- paste(nm1, "*", nm2, "*", geneNm, sep = "")
          id2 <- paste(nm2, "*", nm1, "*", geneNm, sep = "")
          wt <- 1
          weights[[id1]] <- weights[[id2]] <- 1

          lookup[[nm1]][[nm2]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm,
            weight = 1
          )
          lookup[[nm2]][[nm1]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm,
            weight = 1
          )
          # only from large module (src) to small (no need for reverse)
          if (de.num[nm1] >= de.num[nm2]) {
            data[[item.count]] <- list(
              fromModule = nm1,
              toModule = nm2,
              arc = geneNm,
              weight = wt,
              fromColor = as.character(colors[nm1]),
              toColor = as.character(colors[nm2])
            )
          } else {
            data[[item.count]] <- list(
              fromModule = nm2,
              toModule = nm1,
              arc = geneNm,
              weight = wt,
              fromColor = as.character(colors[nm2]),
              toColor = as.character(colors[nm1])
            )
          }
        }
      }
    }
  }

  # need to remove 3rd level list name
  for (m in 1:length(dataNms)) { # need to remove itself
    nm1 <- dataNms[m]
    for (n in m:length(dataNms)) {
      nm2 <- dataNms[n]
      nlist <- lookup[[nm1]][[nm2]]
      lookup[[nm1]][[nm2]] <- lookup[[nm2]][[nm1]] <- unname(nlist)
    }
  }

  chordData <- list(
    arcs = arcs,
    links = links.new,
    weights = weights,
    data = data,
    colors = colors,
    labels = labels
  )

  return(
    list(
      chordData = chordData,
      lookup = lookup
    )
  )
}

PrepareMetaChordData <- function() {
  BHth <- GlobalCutOff$BHth
  logFC <- GlobalCutOff$logFC
  sel.nms <- names(mdata.all)[mdata.all == 1]
  row.nm <- length(sel.nms)

  # update inmex.ind for only selected dataset
  hit.inx <- names(inmex.ind) %in% sel.nms
  inmex.ind <- inmex.ind[hit.inx]
  dataNms <- names(inmex.ind)
  newNms <- substring(dataNms, 0, nchar(dataNms) - 4)

  if (meta.selected) {
    newNms <- c(newNms, "meta")
    row.nm <- row.nm + 1
  }
  inmex.meta <- readRDS("inmex_meta.rds")
  hit.mat <- matrix(0, nrow = row.nm, ncol = nrow(inmex.meta$data))
  nm.mat <- matrix(NA, nrow = row.nm, ncol = nrow(inmex.meta$data))
  shared.ids <- rownames(inmex.meta$data)
  shared.sbls <- inmex.meta$gene.symbls

  colnames(hit.mat) <- colnames(nm.mat) <- shared.sbls
  rownames(hit.mat) <- rownames(nm.mat) <- newNms

  chord.vec.nms <- newNms
  chord.res <- list()
  tot.count <- 0 # need to see how many chords
  for (i in 1:length(inmex.ind)) {
    dataSet <- readRDS(sel.nms[i])
    nm <- newNms[i]
    res.mat <- inmex.ind[[i]][shared.ids, ]
    hit.inx <- abs(res.mat[, 1]) >= logFC & res.mat[, 2] <= dataSet$pval
    hit.mat[i, hit.inx] <- 1
    nm.mat[i, hit.inx] <- nm
    chord.res[[nm]] <- shared.sbls[hit.inx]
    tot.count <- tot.count + sum(hit.inx)
  }

  if (meta.selected) { # add meta
    i <- length(inmex.ind) + 1
    hit.inx <- shared.ids %in% meta.stat$de
    hit.mat[i, hit.inx] <- 1
    nm.mat[i, hit.inx] <- newNms[i]
    chord.res[[newNms[i]]] <- shared.sbls[hit.inx]
    tot.count <- tot.count + sum(hit.inx)
  }

  if (tot.count > 2000) {
    current.msg <<- paste(
      "Chord diagrams is effective to display relationships for less than 1000 genes. The results contain", tot.count,
      "of genes (max. allowed: 2000). You can try Venn diagram or adjust threshold to select most significant genes."
    )
    print(current.msg)
    return(NULL)
  }

  # keep gene with at least 1 hit across all datasets
  keep.inx <- apply(hit.mat, 2, sum) > 0
  hit.mat <- hit.mat[, keep.inx]
  nm.mat <- nm.mat[, keep.inx]
  shared.sbls <- shared.sbls[keep.inx]

  ############## Links#######
  # links & arcs data
  links <- arcs <- list()

  # arcs data
  de.num <- apply(hit.mat, 1, sum)
  de.prct <- as.numeric(de.num / sum(de.num))

  for (i in 1:length(newNms)) {
    arcs[[i]] <- list(
      name = newNms[i],
      label = newNms[i],
      size = as.numeric(de.num[i]),
      val = de.prct[i]
    )
  }

  # now, get unique names by each inserting in the gene names
  lnk.genes <- colnames(nm.mat)
  nm.mat <- apply(
    nm.mat, 1,
    function(x) {
      gd.hits <- !is.na(x)
      x[gd.hits] <- paste(x[gd.hits], "*", lnk.genes[gd.hits], sep = "")
      x
    }
  )

  nm.mat <- unname(t(nm.mat))
  for (l in 1:nrow(nm.mat)) {
    x <- nm.mat[l, ]
    hit.inx <- !is.na(x)
    links$name <- c(links$name, x[hit.inx])

    y.mat <- nm.mat[, hit.inx, drop = F]
    new.mat <- apply(
      y.mat, 2,
      function(d) {
        myd <- d[-l]
        if (all(is.na(myd))) { # link to itself?
          d[l]
        } else {
          myd[!is.na(myd)]
        }
      }
    )
    links$imports <- c(links$imports, new.mat)
  }

  # need to reformat
  links.new <- list()
  for (l in 1:length(links$name)) {
    impt <- links$imports[[l]]
    if (length(impt) == 1) {
      impt <- as.matrix(impt)
    }
    links.new[[l]] <- list(
      name = links$name[l],
      imports = impt
    )
  }

  ###################
  # colors & labels
  library(RColorBrewer)
  if (length(newNms) > 9) {
    colors <- brewer.pal(length(newNms), "Set3")
  } else if (length(newNms) > 2) {
    colors <- brewer.pal(length(newNms), "Dark2")
  } else {
    colors <- c("#7570B3", "#D95F02")
  }
  labels <- rep("", length(newNms))
  names(labels) <- names(colors) <- newNms

  ###################
  # links done!, now get pvalues/data for every pairs of chord
  dataNms <- rownames(hit.mat)
  geneNms <- colnames(hit.mat)
  entrezIDs <- names(shared.sbls)
  weights <- data <- list()

  # also need to setup lookup json
  lookup <- list()
  for (i in 1:length(dataNms)) {
    nlst <- vector(length = length(dataNms) - 1, mode = "list")
    nm <- dataNms[i]
    names(nlst) <- dataNms[-i]
    lookup[[nm]] <- nlst
  }

  item.count <- 0
  for (i in 1:ncol(nm.mat)) {
    geneNm <- geneNms[i]
    entrez <- entrezIDs[i]
    hit.inx <- !is.na(nm.mat[, i])
    datNms <- dataNms[hit.inx]

    if (length(datNms) == 1) {
      item.count <- item.count + 1
      id <- paste(datNms, "*", datNms, "*", geneNm, sep = "")
      weights[[id]] <- 1
      lookup[[datNms]][[datNms]][[geneNm]] <- list(
        entrez = entrez,
        genesymbol = geneNm
      )
      data[[item.count]] <- list(
        fromModule = datNms,
        toModule = datNms,
        arc = geneNm,
        fromColor = as.character(colors[datNms]),
        toColor = as.character(colors[datNms])
      )
    } else {
      for (m in 1:(length(datNms) - 1)) {
        nm1 <- datNms[m]
        for (n in (m + 1):length(datNms)) {
          item.count <- item.count + 1
          nm2 <- datNms[n]
          id1 <- paste(nm1, "*", nm2, "*", geneNm, sep = "")
          id2 <- paste(nm2, "*", nm1, "*", geneNm, sep = "")
          weights[[id1]] <- weights[[id2]] <- 1

          lookup[[nm1]][[nm2]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm
          )
          lookup[[nm2]][[nm1]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm
          )
          # only from large module (src) to small (no need for reverse)
          if (de.num[nm1] >= de.num[nm2]) {
            data[[item.count]] <- list(
              fromModule = nm1,
              toModule = nm2,
              arc = geneNm,
              fromColor = as.character(colors[nm1]),
              toColor = as.character(colors[nm2])
            )
          } else {
            data[[item.count]] <- list(
              fromModule = nm2,
              toModule = nm1,
              arc = geneNm,
              fromColor = as.character(colors[nm2]),
              toColor = as.character(colors[nm1])
            )
          }
        }
      }
    }
  }

  # need to remove 3rd level list name
  for (m in 1:length(dataNms)) { # need to remove itself
    nm1 <- dataNms[m]
    for (n in m:length(dataNms)) {
      nm2 <- dataNms[n]
      nlist <- lookup[[nm1]][[nm2]]
      lookup[[nm1]][[nm2]] <- lookup[[nm2]][[nm1]] <- unname(nlist)
    }
  }

  chordData <- list(
    arcs = arcs,
    links = links.new,
    data = data,
    weights = weights,
    colors = colors,
    labels = labels
  )

  chord.vec.nms <<- chord.vec.nms
  chord.res <<- chord.res

  return(
    list(
      chordData = chordData,
      lookup = lookup
    )
  )
}

PerformChordEnrichment <- function(file.nm, fun.type, IDs) {
  gene.vec <- unlist(strsplit(IDs, "; "))
  sym.vec <- doEntrez2SymbolMapping(gene.vec)
  names(gene.vec) <- sym.vec
  res <- PerformEnrichAnalysis(file.nm, fun.type, gene.vec)
  return(res)
}

PrepareListChordData <- function() {
  all.enIDs <- NULL
  newDat <- list()
  tot.count <- 0
  all.nms <- names(mdata.all)[mdata.all == 1]
  for (i in 1:length(all.nms)) {
    dataNm <- all.nms[i]
    dataSet <- readRDS(dataNm)
    gene.mat <- dataSet$prot.mat

    # convert to entrez
    expr.val <- gene.mat[, 1]
    en.ids <- rownames(gene.mat)

    names(expr.val) <- en.ids
    newDat[[dataNm]] <- expr.val

    all.enIDs <- c(all.enIDs, en.ids)
    tot.count <- tot.count + nrow(gene.mat)

    if (tot.count > 2000) {
      current.msg <<- paste(
        "Chord diagrams is effective to display relationships for less than 1000 items. The results contain", tot.count,
        "of genes (max. allowed: 2000). You can try Venn diagram instead."
      )
      return(NULL)
    }
  }
  PrepareChordDataFromList(newDat, unique(all.enIDs))
}

CalculateDEgeneSet <- function(nms, operation, refNm, filenm) {
  nms <- strsplit(nms, ";")[[1]]
  if (anal.type == "metadata") {
    com.smbls <- PerformSetOperation_Data(nms, operation, refNm)
  } else {
    com.smbls <- PerformSetOperation_List(nms, operation, refNm)
  }

  sink(filenm)
  cat(toJSON(com.smbls))
  sink()
}

PerformSetOperation_List <- function(nms, operation, refNm) {
  all.nms <- names(mdata.all)
  include.inx <- all.nms %in% nms
  my.vec <- all.nms[include.inx]
  if (operation == "diff") { # make sure reference is the first
    inx <- which(my.vec == refNm)
    my.vec <- my.vec[-inx]
  }
  com.ids <- NULL
  ids.list <- list()
  for (i in 1:length(my.vec)) {
    dataSet <- readRDS(my.vec[i])
    if (operation == "diff") {
      ids.list[[i]] <- dataSet$GeneAnotDB[, "gene_id"]
      # com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    } else if (is.null(com.ids)) {
      com.ids <- dataSet$GeneAnotDB[, "gene_id"]
    } else {
      if (operation == "intersect") {
        com.ids <- intersect(com.ids, dataSet$GeneAnotDB[, "gene_id"])
      } else if (operation == "union") {
        com.ids <- union(com.ids, dataSet$GeneAnotDB[, "gene_id"])
      }
    }
  }
  if (operation == "diff") {
    dataSet <- readRDS(refNm)
    ids <- unique(unlist(ids.list))
    com.ids <- setdiff(dataSet$GeneAnotDB[, "gene_id"], ids)
  }

  com.ids <- as.character(com.ids[!is.na(com.ids)]) # make sure it is interpreted as name not index
  com.symbols <- doEntrez2SymbolMapping(com.ids)
  names(com.symbols) <- com.ids

  com.symbols <- com.symbols[!is.null(com.symbols)]
  return(com.symbols)
}

PerformSetOperation_Data <- function(nms, operation, refNm) {
  include.inx <- chord.vec.nms %in% nms
  my.vec <- chord.vec.nms[include.inx]
  my.vec <- paste0(my.vec, ".txt")
  refNm <- paste0(refNm, ".txt")
  if (operation == "diff") { # make sure reference is the first
    inx <- which(my.vec == refNm)
    my.vec <- my.vec[-inx]
  }
  com.ids <- NULL
  ids.list <- list()
  for (nm in my.vec) {
    dataSet <- readRDS(nm)
    if (operation == "diff") {
      ids.list[[nm]] <- rownames(dataSet$sig.mat)
      # com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    } else if (is.null(com.ids)) {
      com.ids <- rownames(dataSet$sig.mat)
    } else {
      if (operation == "intersect") {
        com.ids <- intersect(com.ids, rownames(dataSet$sig.mat))
      } else if (operation == "union") {
        com.ids <- union(com.ids, rownames(dataSet$sig.mat))
      }
    }
  }
  if (operation == "diff") {
    dataSet <- readRDS(refNm)
    ids <- unique(unlist(ids.list))
    com.ids <- setdiff(rownames(dataSet$sig.mat), ids)
  }
  com.symbols <- doEntrez2SymbolMapping(com.ids)
  names(com.symbols) <- com.ids
  return(com.symbols)
}

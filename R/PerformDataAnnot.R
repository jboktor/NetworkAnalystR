#' @title PerformDataAnnot
#' @description read in the data and perform gene ID mapping using built in libraries
#' @param org type of organism
#' @param dataType type of data; either "array" or "count"
#' @param idType INVEX supported ID types
#' @param lvlOpt "NA" to keep original, other values will merge original ID to entrez gene IDs
#' @return return the total matched gene number
#' @details unmapped IDs will be retained as original label (i.e. intergenic regions) in further analysis. matchMin is minimal matched probe (%)
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  PerformDataAnnot("hsa", "array", "hgu95av2", "mean")
#'  }
#' }
#' @rdname PerformDataAnnot
#' @export 
PerformDataAnnot <- function(org, dataType, idType, lvlOpt){
  data.org <<- org;
  SetInitLib(org)
  
  dataSet$type <- dataType;
  dataSet$id.orig <- dataSet$id.current <- idType;
  dataSet$annotated <- F;
  # should not contain duplicates, however sanity check
  data.proc <- readRDS("data.proc.rds");
  dataSet$data.anot <- data.proc;
  
  if (org != 'NA' & idType != 'NA'){
    feature.vec <- rownames(data.proc);
    anot.id <- doAnnotation(feature.vec, idType);
    
    #dataSet$annotation <- anot.id; 
    saveRDS(anot.id, "annotation.rds");
    
    hit.inx <- !is.na(anot.id);
    matched.len <- sum(hit.inx);
    perct <- round(matched.len/length(feature.vec),3)*100;
    thresh <- 0.1 # previous value of 0.25 is causing challenges 
    #for datasets like Ppromelas with low annotation quality
    if (matched.len < length(feature.vec)*thresh){
      current.msg <<- paste('Only ', perct, '% ID were matched. You may want to choose another ID type or use default.', sep=""); 
    } else {
      current.msg <<- paste("ID annotation: ", "Total [", length(anot.id), 
                            "] Matched [", matched.len, "] Unmatched [", sum(!hit.inx),"]", collapse="\n");    
      
      if (lvlOpt != 'NA' | idType == "entrez"){
        # do actual summarization to gene level
        
        matched.entrez <- anot.id[hit.inx];
        data.anot <- data.proc[hit.inx,];
        rownames(data.anot) <- matched.entrez;
        current.msg <<- paste(current.msg, "Data is now transformed to gene-level (Entrez) expression.");
        
        dataSet$data.anot <- RemoveDuplicates(data.anot, lvlOpt, quiet=F);
        dataSet$id.current <- "entrez";
        dataSet$annotated <- T; 
      } else {
        current.msg <<- paste(current.msg, "No gene level summarization was performed.");
      }
    }
  } else { # no conversion will be performed
    feature.vec <- rownames(data.proc);
    anot.id = feature.vec
    perct <- 100;
    hit.inx <- !is.na(anot.id);
    matched.len <- length(feature.vec); # dummies
    minLvl <- 1;
    current.msg <<- paste("No annotation was performed. Make sure organism and gene ID are specified correctly!"); 
  }
  # need to save the ids (mixed gene annotation and original id) 
  # in case, users needs to keep unannotated features
  # this need to be updated to gether with data from now on
  dataSet$data.norm <- dataSet$data.anot;
  dataSet <<- dataSet;
  
  saveRDS(dataSet$data.anot, file="orig.data.anot"); # keep original copy, not in mem
  
  totalCount =  sum(colSums(dataSet$data.anot));
  avgCount = sum(colSums(dataSet$data.anot))/ ncol(dataSet$data.anot);
  minCount = min(colSums(dataSet$data.anot))
  maxCount = max(colSums(dataSet$data.anot))
  
  if(length(dataSet$meta.info)==1){
    lvls = paste(levels(dataSet$meta.info[,1]),collapse="; ")
  }else{
    conc1 = paste0("<b>", colnames(dataSet$meta.info)[1], "</b>", ": ", paste(levels(dataSet$meta.info[,1]), collapse="; "))
    conc2 = paste0("<b>", colnames(dataSet$meta.info)[2], "</b>", ": ", paste(levels(dataSet$meta.info[,2]), collapse="; "))
    lvls = paste("Two factors found -", conc1, conc2)
  }
  summaryVec <<- c(matched.len, perct, length(anot.id), sum(!hit.inx), ncol(dataSet$data.anot), ncol(dataSet$meta.info), sprintf("%4.2e", signif(totalCount ,3)), sprintf("%4.2e",signif(avgCount, 3)), sprintf("%4.2e",signif(minCount, 3)), sprintf("%4.2e",signif(maxCount,3)), lvls)  
  return(matched.len);   
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param id.vec PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doAnnotation
#' @export 
doAnnotation <- function(id.vec, idType){
  if(idType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene","emblprotein", "embltranscript", "orfid", "tair", "wormbase")){
    anot.id <- doGeneIDMapping(id.vec, idType);
  }else{
    anot.id <- doProbeMapping(id.vec, idType);
    names(anot.id) <- id.vec;
  }
  return(anot.id);        
}

# from probe ID to entrez ID 
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param probe.vec PARAM_DESCRIPTION
#' @param platform PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doProbeMapping
#' @export 
doProbeMapping <- function(probe.vec, platform){
  lib.path="data/"
  platform.path <- paste(lib.path,  data.org, "/", platform, ".qs", sep="");
  probe.map <- qs::qread(platform.path);
  if(is.null(probe.vec)){
    entrez <- probe.map[, "entrez"];
  }else{
    hit.inx <- match(probe.vec, probe.map[, "probe"]);
    entrez <- probe.map[hit.inx, "entrez"];
  }
  return(entrez);
}

# mapping between genebank, refseq and entrez
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doGeneIDMapping
#' @export 
doGeneIDMapping <- function(q.vec, type){
  if(is.null(q.vec)){
    db.map <-  queryGeneDB("entrez", data.org);
    q.vec <- db.map[, "gene_id"];
    type = "entrez";
  }
  if(type == "symbol"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    # note, some ID can have version number which is not in the database
    # need to strip it off NM_001402.5 => NM_001402
    q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
    q.vec <- q.mat[,1];
    if(type == "genbank"){
      db.map <-  queryGeneDB("entrez_gb", data.org);
    }else if(type == "refseq"){
      db.map <-  queryGeneDB("entrez_refseq", data.org);
    }else if(type == "emblgene"){
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else if(type == "embltranscript"){
      db.map <-  queryGeneDB("entrez_embl_transcript", data.org);
    }else if(type == "emblprotein"){
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }else if(type == "orfid"){ # only for yeast
      db.map <-  queryGeneDB("entrez_orf", data.org);
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else if(type == "tair"){ # only for ath
      db.map <-  queryGeneDB("tair", data.org);
    }else if(type == "wormbase"){ # only for cel
      db.map <-  queryGeneDB("entrez_wormbase", data.org);
    }else{
      print("Unknown data type");
      return(0);
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
  }
  entrezs=db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  rm(db.map, q.vec); gc();
  return(entrezs);
}

#' @title RemoveDuplicates
#' @description given a data with duplicates, dups is the one with duplicates
#' @param data a data.frame
#' @param lvlOpt "NA" to keep original, other values will merge original ID to entrez gene IDs
#' @param quiet message to user, Default: T
#' @return data.frame without duplicates
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  RemoveDuplicates(data.anot, lvlOpt, quiet=F)
#'  }
#' }
#' @rdname RemoveDuplicates
#' @export 
RemoveDuplicates <- function(data, lvlOpt, quiet=T){
  
  all.nms <- rownames(data);
  colnms <- colnames(data);
  dup.inx <- duplicated(all.nms);
  dim.orig  <- dim(data);
  data <- apply(data, 2, as.numeric); # force to be all numeric
  dim(data) <- dim.orig; # keep dimension (will lost when only one item) 
  rownames(data) <- all.nms;
  colnames(data) <- colnms;
  if(sum(dup.inx) > 0){
    uniq.nms <- all.nms[!dup.inx];
    uniq.data <- data[!dup.inx,,drop=F];
    
    dup.nms <- all.nms[dup.inx];
    uniq.dupnms <- unique(dup.nms);
    uniq.duplen <- length(uniq.dupnms);
    
    for(i in 1:uniq.duplen){
      nm <- uniq.dupnms[i];
      hit.inx.all <- which(all.nms == nm);
      hit.inx.uniq <- which(uniq.nms == nm);
      
      # average the whole sub matrix 
      if(lvlOpt == "mean"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
      }else if(lvlOpt == "median"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
      }else if(lvlOpt == "max"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
      }
    }
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(uniq.data);
  }else{
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(data);
  }
} 

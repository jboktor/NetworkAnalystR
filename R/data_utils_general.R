##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# init resources for analysis
Init.Data<-function(path="../../"){
  .on.public.web <<- FALSE;
  
  partialToBeSaved <<- c("Rload.RData", "Rhistory.R")
  selectedFactorInx <<- 1;
  data.idType <<- "NA";
  chord_count <<-0;
  netUploadU <<-0;
  net.stats <<- as.data.frame(matrix(0, ncol = 3, nrow = 1));
  summaryVec <<- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "NA"); 
  enr.mat <<- NULL;
  numOfLists <<- 1;
  rankOptGlobal <<- "pval";
  data.org <<- "hsa";
  keggpw.count <<- 0;
  pvalu <<- 0.05;

  dataSet <- list(annotated=FALSE);
  dataSet$jsonNms <- list()
  dataSet <<- dataSet
  if(file.exists("/home/glassfish/sqlite/")){
    sqlite.path <<- "/home/glassfish/sqlite/";  #public server
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/"; #xia local
  }else if(file.exists("/home/zzggyy/Downloads/netsqlite/")){
    sqlite.path <<-"/home/zzggyy/Downloads/netsqlite/"; #zgy local
    genesdb.path <<-"/home/zzggyy/Downloads/GeneID_25Species_JE/GeneID_25Species_JE/"
  }else if(file.exists("/home/soufanom/Documents/Projects/Lechang/gene-id-mapping/")){
    sqlite.path <<- "~/Documents/Projects/sqlite/networkanalyst/"; #soufanom local
    genesdb.path <<- "/home/soufanom/Documents/Projects/Lechang/gene-id-mapping/"
  }else if(file.exists("~/Documents/Projects/networkanalyst/")){
    sqlite.path <<- "~/Documents/Projects/networkanalyst/"; #soufanom laptop
    genesdb.path <<- "~/Documents/Projects/gene-id-mapping/"
  }else if(file.exists("/home/le/sqlite/networkanalystdatabase/")){
    sqlite.path <<- "/home/le/sqlite/networkanalystdatabase/"; #le local
    genesdb.path <<- "/home/le/sqlite/GeneID_25Species_JE/"
  }
  
  lib.path <<- paste0(path, "data/");
  data.org <<- NULL;
  module.count <<- 0;
  msg.vec <<- vector(mode="character");
  current.msg <<- "";
  
  # preload some general package
  require('Cairo');
  CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
  require('igraph');
  print("called networkanalyst init!");
}

# genelist, onedata, metadata
# also set up or clear the other global objects
SetAnalType <- function(analType){
  anal.type <<- analType;
  mdata.all <<- list();
  meta.selected <<- TRUE;
  meta.upload <<- FALSE; # when upload merged data from meta-analysis b4
  if(analType == "metadata"){
    partialToBeSaved <<- c(partialToBeSaved, "inmex_meta.rds")
  }
}

SetNetType <- function(netType){
  net.type <<- netType;
}

# When multiple genelists/datasets/results, record their name and save the data as .RDS file
# a) Current dataSet object
# Note, the memory will only contain one dataSet object. By default, the last one will be the current dataSet object;
# Users can switch this (from the interface) to specify which data is load into memory (dataSet object)
# b) Include for certain analysis
# For chord and heatmap analysis, users can do multiple selection (include)
# All datasets are selected by default (1 for selected, 0 for unselected)

# note, dataSet need to have "name" property
RegisterData <- function(dataSet){
  dataName <- dataSet$name;
  saveRDS(dataSet, file=dataName);
  if(!is.null(dataSet$data.raw)){ # save memory for meta-analysis mode
    dataSet$data.raw <- NULL;
  }
  dataSet <<- dataSet;
  if(anal.type == "metadata"){
    mdata.all[[dataName]] <<- 1;
  }else{
    mdata.all <<- lapply(mdata.all, function(x){ x <- 0;});
    mdata.all[[dataName]] <<- 1;
  }
  return(1);
}

# only for switching single expression data results
SetCurrentData <- function(nm){
  if(dataSet$name != nm){
    dataSet <<- readRDS(nm);
  }
  return(1);
}

# remove data object, the current dataSet will be the last one by default
RemoveData <- function(dataName){
  if(!is.null(mdata.all[[dataName]])){
    mdata.all[[dataName]] <<- NULL;
  }
}

# users can select one or more data for analysis
# note, we use 1 to indicate this is selected
# and by default is all selected.
SelectDataSet <- function(){
  if(!exists('nm.vec')){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <- 1;
    }else{
      mdata.all[[nm]] <- 0;
    }
  }
  mdata.all <<- mdata.all;
  
  if(anal.type == "metadata"){
    if("meta_dat" %in% nm.vec){
      meta.selected <<- TRUE;
    }else{
      meta.selected <<- FALSE;
    }
  }
  
  rm('nm.vec', envir = .GlobalEnv);
  return(1);
}

GetAllDataNames <- function(){
  names(mdata.all);
}

SetOrganism <- function(org){
  data.org <<- org;
  init.lib <<- "kegg"
}

SetSelectedFactorInx <- function(inx){
  selectedFactorInx <<- inx;
}

SetSelNetDataset <- function(type){
  selectedNetDataset <<- type;
}

SetRankingMetric <- function(opt){
  rankOptGlobal <<- opt;
}


SetListNms <- function(){
  newDat <- list();
  tot.count <- 0;
  listSizes <- list();
  
  # convert to entrez
  if(anal.type == "metadata"){
    inmex.meta <- readRDS("inmex_meta.rds");
    en.ids <- rownames(inmex.meta$data);
    nm = "meta_data"
  }else{
    en.ids <- rownames(dataSet$resTable)
    nm = "dataSet"
  }
  names(en.ids) <- doEntrez2SymbolMapping(en.ids)
  
  listSizes[[1]] <- list(
    name = nm,
    label = nm,
    size = length(en.ids)
  );
  
  list.genes <<- en.ids;
  listSizes <<- listSizes;
}

SetInitLib <- function(org){
  init.lib <<- "kegg"
}

GetDataListNames <- function(){
  return(names(mdata.all));
}

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformEnrichAnalysis <- function(file.nm, fun.type, ora.vec){
  # prepare lib
  current.geneset <- .loadEnrLib(fun.type)
  
  # prepare query
  ora.nms <- names(ora.vec);
  
  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes
  current.universe <- unique(unlist(current.geneset)); 
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
  
  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.geneset, 
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  
  saveRDS(hits.query, "hits_query.rds");
  
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);
  
  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] <= 0.05;
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }else{
    return(0);
  }
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  enr.mat <<- res.mat
  current.msg <<- "Functional enrichment analysis was completed";
  
  # write json
  require(RJSONIO);
  fun.anot = hits.query; 
  total = resTable[,2]; if(length(total) ==1) { total <- matrix(total) };
  fun.pval = resTable[,5]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj = resTable[,6]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  hit.num = resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(current.setids[names(fun.anot)]); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total
  );
  json.mat <- toJSON(json.res, .na='null');
  json.nm <- paste(file.nm, ".json", sep="");
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  fun.hits <<- hits.query;
  fun.pval <<- resTable[,5];
  hit.num <<- resTable[,4];
  csv.nm <- paste(file.nm, ".csv", sep="");    
  write.csv(resTable, file=csv.nm, row.names=F);
  return(1);
}

GetCurrentJson <-function(type){
  return(dataSet$jsonNms[[type]]);
}

GetFilesToBeSaved <-function(naviString){
  return(unique(partialToBeSaved));
}

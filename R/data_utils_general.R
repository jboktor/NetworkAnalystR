##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# init resources for analysis
Init.Data<-function(path="../../"){
  selectedFactorInx <<- 1;
  fundbu <<- "kegg";
  chord_count <<-0;
  numOfLists <<- 1;
  rankOptGlobal <<- "pval";
  initmetaloading <<- TRUE;
  data.org <<- "hsa";
  keggpw.count <<- 0;
  pvalu <<- 0.05;
  dataSet <<- list(annotated=FALSE); 
  if(file.exists("/home/glassfish/sqlite/networkanalyst/")){
    sqlite.path <<- "/home/glassfish/sqlite/networkanalyst/";  #public server
    # disable parallel prcessing for DESeq2/edgeR
    library(BiocParallel);
    register(SerialParam());
  }else if(file.exists("/Users/xia/Dropbox/sqlite/networkanalyst/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/networkanalyst/"; #xia local
  }else if(file.exists("/home/zzggyy/Downloads/netsqlite/")){
    sqlite.path <<-"/home/zzggyy/Downloads/netsqlite/"; #zgy local
  }else if(file.exists("/home/soufanom/Documents/Projects/sqlite/networkanalyst/")){
    sqlite.path <<- "/home/soufanom/Documents/Projects/sqlite/networkanalyst/"; #soufanom local
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
SelectData <- function(){
  if(!exists('nm.vec')){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <<- 1;
    }else{
      mdata.all[[nm]] <<- 0;
    }
  }
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

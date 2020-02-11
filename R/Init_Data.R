##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
# init resources for analysis
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param path PARAM_DESCRIPTION, Default: '../../'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Init.Data
#' @export 
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
    genesdb.path <<- "/home/glassfish/sqlite/"
    # disable parallel prcessing for DESeq2/edgeR
    library(BiocParallel);
    register(SerialParam());
  }else if(file.exists("/Users/xia/Dropbox/sqlite/networkanalyst/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/networkanalyst/"; #xia local
    genesdb.path <<- "/Users/xia/Dropbox/sqlite/"
  }else if(file.exists("/home/zzggyy/Downloads/netsqlite/")){
    sqlite.path <<-"/home/zzggyy/Downloads/netsqlite/"; #zgy local
    genesdb.path <<-"/home/zzggyy/Downloads/netsqlite/"
  }else if(file.exists("~/Downloads/sqlite")){
    sqlite.path <<- "~/Documents/Projects/sqlite/networkanalyst/"; #soufanom local
    genesdb.path <<- "~/Downloads/sqlite/"
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

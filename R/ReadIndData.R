##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
# read individual data from user, data.type: array
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param data.type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ReadIndData
#' @export 
ReadIndData <- function(dataName, data.type){
  
  current.msg <<- "";
  meta.upload <<- FALSE; # upload data to resume
  dataSet <- .readTabData(dataName);
  
  # now need to remove low quality samples and genes
  data <- dataSet$data;
  meta.info <- dataSet$meta.info;
  
  smpl.num <- ncol(data);
  gene.num <- nrow(data);
  
  # remove smpls/exp with over half missing value
  good.inx<-apply(is.na(data), 2, sum)/nrow(data)<0.6;
  smpl.msg <- "";
  if(sum(!good.inx)>0){
    smpl.msg <- paste(sum(!good.inx), "low quality samples(>60% missing) removed.");
    print(smpl.msg);
    
    data <- data[,good.inx];
    if(ncol(data)/smpl.num < 0.5){
      current.msg <<- paste(smpl.msg, "Low quality data rejected!");;
      return("F");
    }
    
    # update meta information
    meta.info <- meta.info[good.inx, , drop=F];
  }
  
  if(ncol(data) < 4){
    current.msg <<- paste(smpl.msg, "The sample # (", ncol(data), ") is too small.");
    return("F");
  }
  
  # genes with 75% NA will be removed
  gd.inx<-apply(is.na(data), 1, sum)/ncol(data)<0.75;
  feat.msg <- "";
  if(sum(!gd.inx) > 0){
    data <- data[gd.inx,];
    feat.msg <- paste(sum(!gd.inx), "low quality genes (>75% missing) removed");
    if(nrow(data)/gene.num < 0.25){
      current.msg <<- paste(feat.msg, "Low quality data rejected.");
      return("F");
    }
    print(feat.msg);
  }
  
  if(nrow(data) < 10){ 
    current.msg <<- paste(feat.msg, "The gene# (", nrow(data), ") is too small (<10).");
    return("F");
  }
  
  # make an copy, only for testing different normalization
  dataSet$data.raw <- data;
  dataSet$data <- data;
  dataSet$type <- data.type;
  dataSet$meta.info <- meta.info;
  dataName <- dataSet$name;
  res <- RegisterData(dataSet);
  if(res == 1){
    return(dataName);
  }else{
    current.msg <<- paste("Cannot add data: ", dataName, ". ", current.msg, sep="");
    return("F");
  }
}

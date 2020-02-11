#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nms PARAM_DESCRIPTION
#' @param operation PARAM_DESCRIPTION
#' @param refNm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformSetOperation_List
#' @export 
PerformSetOperation_List <- function(nms, operation, refNm){
  all.nms <- names(mdata.all);
  include.inx <- all.nms %in% nms;
  my.vec <- all.nms[include.inx];
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  for(i in 1:length(my.vec)){
    dataSet <- readRDS(my.vec[i]);
    if(operation == "diff"){
      ids.list[[i]]=dataSet$GeneAnotDB[,"gene_id"]
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- dataSet$GeneAnotDB[,"gene_id"];
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }else if(operation == "union"){
        com.ids <- union(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }
    }
  }
  if(operation == "diff"){
    dataSet <- readRDS(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(dataSet$GeneAnotDB[,"gene_id"], ids);
  }
  
  com.ids <- as.character(com.ids[!is.na(com.ids)]); # make sure it is interpreted as name not index
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) = com.ids
  
  com.symbols<-com.symbols[!is.null(com.symbols)];
  return(com.symbols);
}

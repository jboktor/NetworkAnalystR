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
#' @rdname PerformSetOperation_Data
#' @export 
PerformSetOperation_Data <- function(nms, operation, refNm){
  include.inx <- chord.vec.nms %in% nms;
  my.vec <- chord.vec.nms[include.inx];
  my.vec <- paste0(my.vec,".txt")
  refNm = paste0(refNm, ".txt")
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  for(nm in my.vec){
    dataSet = readRDS(nm);
    if(operation == "diff"){
      ids.list[[nm]]=rownames(dataSet$sig.mat);
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- rownames(dataSet$sig.mat);
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, rownames(dataSet$sig.mat));
      }else if(operation=="union"){
        com.ids <- union(com.ids, rownames(dataSet$sig.mat));
      }
    }
  }
  if(operation == "diff"){
    dataSet <- readRDS(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(rownames(dataSet$sig.mat), ids);
  } 
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) = com.ids;
  return(com.symbols);
}

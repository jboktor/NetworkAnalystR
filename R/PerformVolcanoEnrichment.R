##################################################
## R script for NetworkAnalyst
## Description: Functions for heatmaps
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param IDs PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @param inx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformVolcanoEnrichment
#' @export 
PerformVolcanoEnrichment<-function(file.nm, fun.type, IDs, type, inx){
  inx = as.numeric(inx)
  if(anal.type == "onedata"){
    if(dataSet$type == "array"){
      sigmat = dataSet$sig.mat
    } else {
      sigmat = dataSet$sig.mat
    }
  }else{
    sigmat = inmex.ind[selDataNm][[1]][which(inmex.ind[selDataNm][[1]][,'Pval'] < as.numeric(pvalu)),];
  }
  sigm <<- sigmat
  
  if(type == "focus"){
    gene.vec <- unlist(strsplit(IDs, "; "));
  }else if(type == "all"){
    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
    gene.vec <- c(gene.vecup, gene.vecdown);
  }else if(type == "up"){
    gene.vec <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
  }else{
    gene.vec <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
  }
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, gene.vec);
  return(res);
}

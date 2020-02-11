#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param IDs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformHeatmapEnrichment
#' @export 
PerformHeatmapEnrichment <- function(file.nm, fun.type, IDs){
  if(IDs=="NA"){
    if(anal.type=="onedata"){
      gene.vec <- rownames(dataSet$sig.mat);
    }else if(anal.type=="metadata"){
      gene.vec <- rownames(meta.mat);
    }else{
      gene.vec <- rownames(all.ent.mat);
    }
  }else{
    gene.vec <- unlist(strsplit(IDs, "; "));
  }
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, gene.vec);
  return(res);
}

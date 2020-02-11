#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param netNm PARAM_DESCRIPTION
#' @param IDs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformListEnrichmentView
#' @export 
PerformListEnrichmentView <- function(file.nm, fun.type, netNm, IDs){
  gene.vec <- unlist(strsplit(IDs, "; "));
  gene.vec <- unique(gene.vec);
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  list.genes <<- gene.vec
  res <- PerformEnrichAnalysis(file.nm, fun.type, list.genes);
  PrepareEnrichNet(netNm, 'list', "mixed");
  return(res);
}

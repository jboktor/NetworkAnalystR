#areas is allname concated
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param areas PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetVennGeneNames
#' @export 
GetVennGeneNames <- function(areas){
  nms <- strsplit(areas, "\\|\\|")[[1]];
  gene.vec <- NULL;
  for(nm in nms){
    gene.vec <- c(gene.vec, vennData[[nm]]);
  }
  gene.vec <- unique(gene.vec);
  # from entrez to symbols
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  venn.genes <<- gene.vec;
  return(paste(unique(sym.vec), collapse="||"));
}

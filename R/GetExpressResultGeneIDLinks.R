#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetExpressResultGeneIDLinks
#' @export 
GetExpressResultGeneIDLinks <- function(){
  ids <- rownames(dataSet$resTable);
  symbs <- doEntrez2SymbolMapping(ids);
  annots <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids, "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformVennEnrichment
#' @export 
PerformVennEnrichment <- function(file.nm, fun.type){
  res <- PerformEnrichAnalysis(file.nm, fun.type, venn.genes);
  return(res);
}

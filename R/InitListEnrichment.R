#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname InitListEnrichment
#' @export 
InitListEnrichment <- function(type){
  GetListEnrGeneNumber();
  res <- PerformEnrichAnalysis(paste0("enrichment_", type), type, list.genes);
  PrepareEnrichNet(paste0('enrichNet_', type), 'list', "mixed");
  return(res)
}

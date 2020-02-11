#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param org PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetOrganism
#' @export 
SetOrganism <- function(org){
  data.org <<- org;
  init.lib <<- "kegg"
}

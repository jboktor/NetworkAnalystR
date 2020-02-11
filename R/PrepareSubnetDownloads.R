#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareSubnetDownloads
#' @export 
PrepareSubnetDownloads <- function(nm){
  g <- ppi.comps[[nm]];
  # need to update graph so that id is compound names rather than ID
  V(g)$name <- as.character(doID2LabelMapping(V(g)$name));
  saveNetworkInSIF(g, nm);
}

# also save to GraphML
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ExportNetwork
#' @export 
ExportNetwork <- function(fileName){
  current.net <- ppi.comps[[current.net.nm]];
  write.graph(current.net, file=fileName, format="graphml");
}

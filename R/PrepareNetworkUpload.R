#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.nm PARAM_DESCRIPTION
#' @param json.nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareNetworkUpload
#' @export 
PrepareNetworkUpload <- function(net.nm, json.nm){
  
  my.ppi <- ppi.comps[[net.nm]];
  nd.nms <- V(my.ppi)$name;
  
  convertIgraph2JSONFromFile(net.nm, json.nm);
  current.net.nm <<- net.nm;
  return(1);
}

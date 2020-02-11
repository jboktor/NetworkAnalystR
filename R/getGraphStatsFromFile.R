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
#' @rdname getGraphStatsFromFile
#' @export 
getGraphStatsFromFile <- function(){
  g <- ppi.comps[[net.nmu]];
  nms <- V(g)$name;
  edge.mat <- get.edgelist(g);
  return(c(length(nms), nrow(edge.mat)));        
}

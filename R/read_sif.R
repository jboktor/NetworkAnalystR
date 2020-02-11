#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sif.file PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION, Default: 'graphNEL'
#' @param directed PARAM_DESCRIPTION, Default: FALSE
#' @param header PARAM_DESCRIPTION, Default: TRUE
#' @param sep PARAM_DESCRIPTION, Default: '	'
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname read.sif
#' @export 
read.sif <- function (sif.file, format = "graphNEL", directed = FALSE, header = TRUE, sep = "\t", ...) {
  
  net <- read.csv(file = sif.file, sep = sep, colClasses = "character", header = header, ...)
  
  # Assume form: node1 linktype node2 side.info..
  if ( ncol(net) > 2 ) { 
    
    # remove NA nodes 
    nas <- apply(net, 1, function (x) {any(is.na(x[c(1,3)]))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }
    
    net <- graph.edgelist(as.matrix(net[, -2]), directed = directed)
    
  } else if ( ncol(net) == 2 ) { # assume form: node1 node2
    
    # remove NA nodes 
    nas <- apply(net, 1, function (x) {any(is.na(x))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }
    
    net <- graph.edgelist(cbind(net[,1],net[,2]), directed = directed)
  }
  
  if (format == "graphNEL") { net <- igraph.to.graphNEL(net) }
  
  return(net);
}

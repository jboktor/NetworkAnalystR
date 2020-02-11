#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param algo PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname UpdateNetworkLayout
#' @export 
UpdateNetworkLayout <- function(algo, filenm){
  current.net <- ppi.comps[[current.net.nm]];
  #pos.xy <- PerformLayOut_mem(current.net.nm, algo);
  pos.xy <- PerformLayOut(current.net.nm, algo);
  nms <- V(current.net)$name;
  nodes <- vector(mode="list");
  for(i in 1:length(nms)){
    nodes[[i]] <- list(
      id=nms[i], 
      x=pos.xy[i,1], 
      y=pos.xy[i,2]
    );
  }
  # now only save the node pos to json
  require(RJSONIO);
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(filenm);
}

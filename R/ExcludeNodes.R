# exclude nodes in current.net (networkview)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nodeids PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ExcludeNodes
#' @export 
ExcludeNodes <- function(nodeids, filenm){
  
  nodes2rm <- strsplit(nodeids, ";")[[1]];
  current.net <- ppi.comps[[current.net.nm]];
  current.net <- delete.vertices(current.net, nodes2rm);
  
  # need to remove all orphan nodes
  bad.vs<-V(current.net)$name[degree(current.net) == 0];
  current.net <- delete.vertices(current.net, bad.vs);
  
  # return all those nodes that are removed 
  nds2rm <- paste(c(bad.vs, nodes2rm), collapse="||");
  
  # update topo measures
  node.btw <- as.numeric(betweenness(current.net));
  node.dgr <- as.numeric(degree(current.net));
  node.exp <- as.numeric(get.vertex.attribute(current.net, name="abundance", index = V(current.net)));
  nms <- V(current.net)$name;
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  lbls <- ppi.net$node.data[hit.inx,2];
  
  nodes <- vector(mode="list");
  for(i in 1:length(nms)){
    nodes[[i]] <- list(
      id=nms[i], 
      label=lbls[i],
      degree=node.dgr[i], 
      between=node.btw[i],
      expr = node.exp[i]
    );
  }
  # now only save the node pos to json
  require(RJSONIO);
  netData <- list(deletes=nds2rm,nodes=nodes);
  sink(filenm);
  cat(toJSON(netData));
  sink();
  
  ppi.comps[[current.net.nm]] <<- current.net;
  UpdateSubnetStats();
  
  # remember to forget the cached layout, and restart caching, as this is now different object (with the same name)
  #forget(PerformLayOut_mem);
  return(filenm);
}

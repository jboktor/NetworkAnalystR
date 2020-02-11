#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nodeids PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ExtractModule
#' @export 
ExtractModule<- function(nodeids){
  set.seed(8574);
  nodes <- strsplit(nodeids, ";")[[1]];
  
  g <- ppi.comps[[current.net.nm]];
  # try to see if the nodes themselves are already connected
  hit.inx <- V(g)$name %in% nodes; 
  gObj <- induced.subgraph(g, V(g)$name[hit.inx]);
  
  # now find connected components
  comps <-decompose.graph(gObj, min.vertices=1);
  
  if(length(comps) == 1){ # nodes are all connected
    g <- comps[[1]];
  }else{
    # extract modules
    paths.list <-list();
    sd.len <- length(nodes);
    for(pos in 1:sd.len){
      paths.list[[pos]] <- get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- V(g)$name[-nds.inxs];
    g <- simplify(delete.vertices(g, nodes2rm));
  }
  nodeList <- get.data.frame(g, "vertices");
  if(nrow(nodeList) < 3){
    return ("NA");
  }
  
  module.count <- module.count + 1;
  module.nm <- paste("module", module.count, sep="");
  colnames(nodeList) <- c("Id", "Label");
  ndFileNm = paste(module.nm, "_node_list.csv", sep="");
  write.csv(nodeList, file=ndFileNm, row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  edgFileNm = paste(module.nm, "_edge_list.csv", sep="");
  write.csv(edgeList, file=edgFileNm, row.names=F, quote=F);
  
  filenm <- paste(module.nm, ".json", sep="");
  
  # record the module 
  ppi.comps[[module.nm]] <<- g;
  UpdateSubnetStats();
  
  module.count <<- module.count
  
  convertIgraph2JSON(module.nm, filenm);
  return (filenm);
}

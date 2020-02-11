# exclude nodes in overall net (network builder)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nodeids PARAM_DESCRIPTION
#' @param id.type PARAM_DESCRIPTION
#' @param vismode PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ExcludeNodesOverall
#' @export 
ExcludeNodesOverall <- function(nodeids, id.type, vismode){
  # all convert to uniprot ID 
  lines <- strsplit(nodeids, "\r|\n|\r\n")[[1]];
  lines<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
  if(vismode != "network"){
    prot.anots <- convertIdToEntrez(lines, id.type);
    nodes2rm <- unique(prot.anots$accession);
  }else{
    prot.anots <- lines
    nodes2rm <- unique(lines);
  }
  
  # now find the overlap
  nodes2rm <- nodes2rm[nodes2rm %in% V(overall.graph)$name];
  g <- delete.vertices(overall.graph, nodes2rm);
  
  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  
  write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);
  
  substats <- DecomposeGraph(g);
  if(!is.null(substats)){
    overall.graph <<- g;
    #forget(PerformLayOut_mem);
    return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

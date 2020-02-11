#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nd.type PARAM_DESCRIPTION
#' @param min.dgr PARAM_DESCRIPTION
#' @param min.btw PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FilterBipartiNet
#' @export 
FilterBipartiNet <- function(nd.type, min.dgr, min.btw){
  
  all.nms <- V(overall.graph)$name;
  edge.mat <- get.edgelist(overall.graph);
  dgrs <- degree(overall.graph);
  nodes2rm.dgr <- nodes2rm.btw <- NULL;
  
  if(nd.type == "gene"){
    hit.inx <- all.nms %in% edge.mat[,1];
  }else if(nd.type=="other"){
    hit.inx <- all.nms %in% edge.mat[,2];
  }else{ # all
    hit.inx <- rep(TRUE, length(all.nms));
  }
  
  if(min.dgr > 0){
    rm.inx <- dgrs <= min.dgr & hit.inx;
    nodes2rm.dgr <- V(overall.graph)$name[rm.inx];
  }
  if(min.btw > 0){
    btws <- betweenness(overall.graph);
    rm.inx <- btws <= min.btw & hit.inx;
    nodes2rm.btw <- V(overall.graph)$name[rm.inx];
  }
  
  nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
  overall.graph <- simplify(delete.vertices(overall.graph, nodes2rm));
  current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
  substats <- DecomposeGraph(overall.graph);
  if(!is.null(substats)){
    overall.graph <<- overall.graph;
    return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

# zero-order network - create ppi nets from only input (seeds)
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
#' @seealso 
#'  \code{\link[igraph]{V}},\code{\link[igraph]{c("delete_vertices", "delete_vertices")}},\code{\link[igraph]{as_data_frame}}
#'  \code{\link[DelayedArray]{simplify}}
#'  \code{\link[AnnotationDbi]{Bimap-toTable}}
#'  \code{\link[S4Vectors]{Vector-class}}
#'  \code{\link[readr]{write_delim}}
#' @rdname BuildSeedProteinNet
#' @export 
#' @importFrom igraph V delete.vertices get.data.frame
#' @importFrom DelayedArray simplify
#' @importFrom AnnotationDbi colnames
#' @importFrom S4Vectors c
#' @importFrom readr write.csv
BuildSeedProteinNet <- function(){
  
  nodes = igraph::V(overall.graph)$name;
  hit.inx <-  nodes %in% seed.proteins;
  nodes2rm <- nodes[!hit.inx];
  g <- DelayedArray::simplify(igraph::delete.vertices(overall.graph, nodes2rm));
  
  nodeList <- igraph::get.data.frame(g, "vertices");
  nodeList <- nodeList[,1:2];
  AnnotationDbi::colnames(nodeList) <- S4Vectors::c("Id", "Label");
  
  readr::write.csv(nodeList, file="orig_node_list.csv", quote=F);
  nd.inx <- ppi.net$node.data[,1] %in% nodeList[,1];
  
  edgeList <- igraph::get.data.frame(g, "edges");
  edgeList <- edgeList[,1:2];
  AnnotationDbi::colnames(edgeList) <- S4Vectors::c("Source", "Target");
  readr::write.csv(edgeList, file="orig_edge_list.csv",  quote=F);
  
  # update ppi.net 
  ppi.net$order = 0;
  ppi.net$node.data <- nodeList;
  ppi.net$edge.data <- edgeList;
  ppi.net <<- ppi.net;
}

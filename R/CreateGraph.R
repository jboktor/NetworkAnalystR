# create igraph from the edgelist saved from graph DB
# and decompose into subnets
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
#' @rdname CreateGraph
#' @export 
CreateGraph <- function(){
  
  require('igraph');
  node.list <- ppi.net$node.data;
  edge.list <- ppi.net$edge.data;
  
  seed.proteins <- ppi.net$seeds;
  overall.graph <- simplify(graph.data.frame(edge.list, directed=FALSE, vertices=node.list));
  
  # add node expression value
  if(ppi.net$db.type == "ppi"){
    ## Temp fix seed.expr all in entrez? 
    ## all converted to new IDs based on PPI from the given organism 
    newIDs <- doPpiIDMapping(names(seed.expr))$accession;    
  }else{# all entrez in mirNet
    newIDs <- names(seed.expr);
  }
  match.index <- match(V(overall.graph)$name, newIDs);
  expr.vals <- seed.expr[match.index];
  #   expr.vec <- abs(expr.vals) # logFC can be negative!!
  expr.vec <- expr.vals;
  names(expr.vec)= V(overall.graph)$name;
  expr.vec <<- expr.vec[!is.na(expr.vec)]
  overall.graph <- set.vertex.attribute(overall.graph, "abundance", index = V(overall.graph), value = expr.vals);
  
  hit.inx <- seed.proteins %in% node.list[,1];
  seed.proteins <<- seed.proteins[hit.inx];
  
  substats <- DecomposeGraph(overall.graph);
  overall.graph <<- overall.graph;
  if(!is.null(substats)){
    return(c(length(seed.genes), length(seed.proteins), nrow(node.list), nrow(edge.list), length(ppi.comps), substats));        
  }else{
    return(0);
  }
}

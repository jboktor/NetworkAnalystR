# for a given graph, obtain the smallest subgraphs that contain
# all the seed nodes. This is acheived by iteratively remove 
# the marginal nodes (degree = 1) that are not in the seeds
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param max.len PARAM_DESCRIPTION, Default: 200
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMinConnectedGraphs
#' @export 
GetMinConnectedGraphs <- function(max.len = 200){
  set.seed(8574);
  # first get shortest paths for all pair-wise seeds
  my.seeds <- seed.proteins;
  sd.len <- length(my.seeds);
  paths.list <-list();
  
  # first trim overall.graph to remove no-seed nodes of degree 1
  dgrs <- degree(overall.graph);
  keep.inx <- dgrs > 1 | (names(dgrs) %in% my.seeds);
  nodes2rm <- V(overall.graph)$name[!keep.inx];
  overall.graph <-  simplify(delete.vertices(overall.graph, nodes2rm));
  
  # need to restrict the operation b/c get.shortest.paths is very time consuming
  # for top max.len highest degrees
  if(sd.len > max.len){
    hit.inx <- names(dgrs) %in% my.seeds;
    sd.dgrs <- dgrs[hit.inx];
    sd.dgrs <- rev(sort(sd.dgrs));
    # need to synchronize all (seed.proteins) and top seeds (my.seeds)
    seed.proteins <- names(sd.dgrs);    
    my.seeds <- seed.proteins[1:max.len];
    sd.len <- max.len;
    current.msg <<- paste("The minimum connected network was computed using the top", sd.len, "seed proteins in the network based on their degrees.");
  }else{
    current.msg <<- paste("The minimum connected network was computed using all seed proteins in the network.");
  }
  # now calculate the shortest paths for 
  # each seed vs. all other seeds (note, to remove pairs already calculated previously)
  for(pos in 1:sd.len){
    paths.list[[pos]] <- get.shortest.paths(overall.graph, my.seeds[pos], seed.proteins[-(1:pos)])$vpath;
  }
  nds.inxs <- unique(unlist(paths.list));
  nodes2rm <- V(overall.graph)$name[-nds.inxs];
  g <- simplify(delete.vertices(overall.graph, nodes2rm));
  
  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);
  
  path.list <- NULL;
  substats <- DecomposeGraph(g);
  if(!is.null(substats)){
    overall.graph <<- g;
    return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

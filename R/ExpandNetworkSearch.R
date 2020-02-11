# 2nd order network, only for PPI
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
#' @rdname ExpandNetworkSearch
#' @export 
ExpandNetworkSearch <- function(){
  if(ppi.net$order == 0){
    SearchNetDB(ppi.net$db.type, ppi.net$table.nm);
    CreateGraph();
  }
  protein.vec <- V(overall.graph)$name;
  table.nm <- ppi.net$table.nm;
  if(db.typeu == "ppi"){
    res <- QueryPpiSQLite(table.nm, protein.vec, ppi.net$require.exp, ppi.net$min.score);
  }else if (db.typeu == "tissuecoex"){
    res <- QueryTissueCoexSQLite(protein.vec);
  }else if (db.typeu == "tissueppi"){
    res <- QueryDiffNetSQLite(protein.vec);
  }else if (db.typeu == "cellcoex"){
    res <- QueryCellCoexSQLite(protein.vec);
  }
  edge.res <- data.frame(Source=res[,1],Target=res[,2]);
  #row.names(edge.res) <- res[,5];
  
  node.ids <- c(res[,1], res[,2])
  node.nms <- c(res[,3], res[,4]);
  node.res <- data.frame(Id=node.ids, Label=node.nms);
  node.res <- node.res[!duplicated(node.res$Id),];
  if(nrow(node.res) < 10000){ # only overwrite if within the range
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    write.csv(node.res, file="orig_node_list.csv", row.names=FALSE);
    ppi.net <<- list(db.type=db.typeu, order=2, seeds=protein.vec, table.nm=table.nm, node.data = node.res, edge.data = edge.res);
  }
  return(nrow(node.res));
}

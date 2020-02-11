# use microservice
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
#' @rdname GetPCSFNet
#' @export 
GetPCSFNet <- function(){
  
  print("Peforming PCSF ....");
  
  library(RSclient);
  rsc <- RS.connect();
  RS.assign(rsc, "my.dir", getwd()); 
  RS.eval(rsc, setwd(my.dir));
  
  dat.out <- list(data=overall.graph, terminals = expr.vec);
  RS.assign(rsc, "dat.in", dat.out); 
  my.fun <- function(){
    require('PCSF');
    edg <- get.edgelist(dat.in$data);
    edg <- as.data.frame(edg);
    edg$V3 <- rep(1, nrow(edg));
    colnames(edg) <- c("from", "to", "cost");
    ppi <- construct_interactome(edg);
    g <- PCSF(ppi, dat.in$terminals, w = 5, b = 100, mu = 0.0005);
    return(g);
  }
  
  RS.assign(rsc, my.fun);
  g <-  RS.eval(rsc, my.fun());
  RS.close(rsc);
  
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

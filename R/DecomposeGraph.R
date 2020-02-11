##################################################
## R script for NetworkAnalyst
## Description: General graph manipulation functions 
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gObj PARAM_DESCRIPTION
#' @param minNodeNum PARAM_DESCRIPTION, Default: 3
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname DecomposeGraph
#' @export 
DecomposeGraph <- function(gObj, minNodeNum = 3){
  # now decompose to individual connected subnetworks
  comps <-decompose.graph(gObj, min.vertices=minNodeNum);
  if(length(comps) == 0){
    current.msg <<- paste("No subnetwork was identified with at least", minNodeNum, "nodes!");
    return(NULL);
  }
  
  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("subnetwork", 1:length(comps), sep="");
  
  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];
  
  # now record
  ppi.comps <<- comps;
  net.stats <<- net.stats;
  
  sub.stats <- unlist(lapply(comps, vcount)); 
  return(sub.stats);
}

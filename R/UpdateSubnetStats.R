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
#' @rdname UpdateSubnetStats
#' @export 
UpdateSubnetStats <- function(){
  old.nms <- names(ppi.comps);
  net.stats <- ComputeSubnetStats(ppi.comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  rownames(net.stats) <- old.nms[ord.inx];
  net.stats <<- net.stats;
}

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
#' @rdname GetNodeBetweenness
#' @export 
GetNodeBetweenness <- function(){
  round(betweenness(overall.graph, directed=F, normalized=F), 2);
}

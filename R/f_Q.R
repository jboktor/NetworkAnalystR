##################################
# functions for estimating Cochran’s Q
##################################
#computes Cochran’s Q gene by gene
#dadj and varadj must be matrices, in which every study is a column,
#every row a gene
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dadj PARAM_DESCRIPTION
#' @param varadj PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname f.Q
#' @export 
f.Q <- function(dadj, varadj){
  w<-1/varadj
  tmp1<-w*dadj
  mu<-rowSums(tmp1)/rowSums(w)
  Q<-rowSums(w*(dadj - mu)^2)
}

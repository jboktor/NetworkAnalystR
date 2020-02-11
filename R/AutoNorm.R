# normalize to zero mean and unit variance
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[BiocGenerics]{mean}},\code{\link[BiocGenerics]{var}}
#' @rdname AutoNorm
#' @export 
#' @importFrom BiocGenerics mean sd
AutoNorm<-function(x){
  (x - BiocGenerics::mean(x))/BiocGenerics::sd(x, na.rm=T);
}

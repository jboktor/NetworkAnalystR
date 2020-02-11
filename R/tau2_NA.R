#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param Q PARAM_DESCRIPTION
#' @param num.studies PARAM_DESCRIPTION
#' @param my.weights PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname tau2.NA
#' @export 
tau2.NA <- function(Q, num.studies, my.weights) {
  vwts <- rowSums(my.weights, na.rm = TRUE)
  tmp2 <- rowSums(my.weights^2, na.rm = TRUE)
  tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
  return(tau2)
}

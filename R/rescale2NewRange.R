# new range [a, b]
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param qvec PARAM_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname rescale2NewRange
#' @export 
rescale2NewRange <- function(qvec, a, b){
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}

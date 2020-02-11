#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param z PARAM_DESCRIPTION
#' @param D PARAM_DESCRIPTION
#' @param d12 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname est.hyper
#' @export 
est.hyper <- function (z, D, d12) 
{
  f2 <- function(d0, D) {
    var(z, na.rm = TRUE) - trigamma(D/2) - trigamma(d0/2)
  }
  lim = f2(100, D)
  if (lim < 0) 
    d0.est <- 100
  if (lim > 0) 
    d0.est <- uniroot(f2, c(1, 100), D = D, extendInt = "yes")$root
  s2.est <- exp(mean(z, na.rm = TRUE) - digamma(D/2) + digamma(d0.est/2) - 
                  log(d0.est/D))
  return(list(d0 = d0.est, s2 = s2.est))
}

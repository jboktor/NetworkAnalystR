#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param d0 PARAM_DESCRIPTION
#' @param D PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname f2
#' @export 
  f2 <- function(d0, D) {
    var(z, na.rm = TRUE) - trigamma(D/2) - trigamma(d0/2)
  }

# overwrite ave, => na.rm=T
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname myave
#' @export 
myave <- function (x, ...) {
  n <- length(list(...))
  if (n) {
    g <- interaction(...)
    split(x, g) <- lapply(split(x, g), mean, na.rm=T)
  }
  else x[] <- FUN(x, na.rm=T)
  x
}

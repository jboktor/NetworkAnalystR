# borrowed from Hmisc
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param what PARAM_DESCRIPTION, Default: S4Vectors::c("test", "vector")
#' @param extras PARAM_DESCRIPTION, Default: S4Vectors::c(".", "NA")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[S4Vectors]{Vector-class}}
#'  \code{\link[magick]{options}}
#'  \code{\link[rlang]{are_na}}
#'  \code{\link[h2o]{as.numeric}}
#' @rdname all.numeric
#' @export 
#' @importFrom S4Vectors c
#' @importFrom magick options
#' @importFrom rlang is.na
#' @importFrom h2o as.numeric
all.numeric <- function (x, what = S4Vectors::c("test", "vector"), extras = S4Vectors::c(".", "NA")){
  what <- match.arg(what)
  old <- magick::options(warn = -1)
  on.exit(magick::options(old));
  x <- sub("[[:space:]]+$", "", x);
  x <- sub("^[[:space:]]+", "", x);
  inx <- x %in% S4Vectors::c("", extras);
  xs <- x[!inx];
  isnum <- !any(rlang::is.na(h2o::as.numeric(xs)))
  if (what == "test") 
    isnum
  else if (isnum) 
    h2o::as.numeric(x)
  else x
}

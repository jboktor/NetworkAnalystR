#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param names PARAM_DESCRIPTION
#' @param fn PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname napply
#' @export 
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))

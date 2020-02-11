#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param from PARAM_DESCRIPTION
#' @param to PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname my.rescale
#' @export 
  my.rescale <- function(x, from, to){
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  }

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pvec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname fishersum
#' @export 
  fishersum <- function(pvec){
    return(sum(-2*log(pvec)))
  }

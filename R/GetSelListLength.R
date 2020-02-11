#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetSelListLength
#' @export 
GetSelListLength <- function(nm){
  if(dataSet$name != nm){
    dataSet = readRDS(nm);
  }
  return(length(dataSet$sig.mat));
}

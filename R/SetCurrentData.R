# only for switching single expression data results
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
#' @rdname SetCurrentData
#' @export 
SetCurrentData <- function(nm){
  if(dataSet$name != nm){
    dataSet <<- readRDS(nm);
  }
  return(1);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetDataDims
#' @export 
GetDataDims <- function(dataName){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  dm <- dim(dataSet$data);
  naNum <- sum(is.na(dataSet$data));
  return(c(dm, naNum));
} 

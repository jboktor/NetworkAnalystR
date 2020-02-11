# remove data object, the current dataSet will be the last one by default 
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
#' @rdname RemoveData
#' @export 
RemoveData <- function(dataName){
  if(!is.null(mdata.all[[dataName]])){
    mdata.all[[dataName]] <<- NULL;
  }
}

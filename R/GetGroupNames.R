# here should first try to load the original data
# the data in the memory could be changed
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
#' @rdname GetGroupNames
#' @export 
GetGroupNames <- function(dataName){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  levels(dataSet$cls);
}

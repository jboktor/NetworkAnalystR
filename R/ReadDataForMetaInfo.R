# read the uploaded data into memory
# return the meta-data information (multiple groups)
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
#' @rdname ReadDataForMetaInfo
#' @export 
ReadDataForMetaInfo<-function(dataName){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  return(colnames(dataSet$meta.info));
}

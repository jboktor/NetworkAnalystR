#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMetaDatasets
#' @export 
GetMetaDatasets<- function(){
  sel.nms <- names(mdata.all)[mdata.all==1];
  return(sel.nms);
}

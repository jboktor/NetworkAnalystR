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
#' @rdname GetMetaResultColNames
#' @export 
GetMetaResultColNames<-function(){
  sel.nms <- names(mdata.all)[mdata.all==1];
  c(substring(sel.nms, 0, nchar(sel.nms)-4), colnames(meta.mat));
}

# genelist, onedata, metadata
# also set up or clear the other global objects
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param analType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetAnalType
#' @export 
SetAnalType <- function(analType){
  anal.type <<- analType;
  mdata.all <<- list(); 
  meta.selected <<- TRUE;
  meta.upload <<- FALSE; # when upload merged data from meta-analysis b4
}

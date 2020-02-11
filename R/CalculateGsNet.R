#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param netNm PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @param mType PARAM_DESCRIPTION
#' @param db PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname CalculateGsNet
#' @export 
CalculateGsNet <- function(name, netNm, type, mType, db){
  res = PerformMetaPathCombine(name, netNm, "pval", db, mType,0.05)
  return(1)
}

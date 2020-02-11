#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nms PARAM_DESCRIPTION
#' @param operation PARAM_DESCRIPTION
#' @param refNm PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[data.table]{tstrsplit}}
#'  \code{\link[R.utils]{Non-documented objects}}
#'  \code{\link[jsonlite]{toJSON, fromJSON}}
#' @rdname CalculateDEgeneSetEnr
#' @export 
#' @importFrom data.table strsplit
#' @importFrom R.utils cat
#' @importFrom jsonlite toJSON
CalculateDEgeneSetEnr <- function(nms, operation, refNm, filenm){
  nms <- data.table::strsplit(nms, ";")[[1]];
  if(anal.type == "metadata" || anal.type == "onedata"){
    com.smbls <- PerformSetOperation_DataEnr(nms, operation, refNm);
  }else{
    com.smbls <- PerformSetOperation_ListEnr(nms, operation, refNm);
  }
  
  sink(filenm);
  R.utils::cat(jsonlite::toJSON(com.smbls));
  sink();
}

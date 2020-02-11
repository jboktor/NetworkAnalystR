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
#' @rdname CalculateDEgeneSet
#' @export 
#' @importFrom data.table strsplit
#' @importFrom R.utils cat
#' @importFrom jsonlite toJSON
CalculateDEgeneSet <- function(nms, operation, refNm, filenm){
  nms <- data.table::strsplit(nms, ";")[[1]];
  if(anal.type == "metadata"){
    com.smbls <- PerformSetOperation_Data(nms, operation, refNm);
  }else{
    com.smbls <- PerformSetOperation_List(nms, operation, refNm);
  }
  
  sink(filenm);
  R.utils::cat(jsonlite::toJSON(com.smbls));
  sink();
}

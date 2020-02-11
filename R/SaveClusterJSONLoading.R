#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileNm PARAM_DESCRIPTION
#' @param clustOpt PARAM_DESCRIPTION
#' @param nb PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SaveClusterJSONLoading
#' @export 
SaveClusterJSONLoading <- function(fileNm, clustOpt, nb){
  if(anal.type == "onedata"){
    SaveExpressClusterLoadingJSON(fileNm, clustOpt, nb);
  }else{
    SaveMetaClusterLoadingJSON(fileNm, clustOpt, nb);
  }
}

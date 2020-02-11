#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileNm PARAM_DESCRIPTION
#' @param clustOpt PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SaveClusterJSON
#' @export 
SaveClusterJSON <- function(fileNm, clustOpt, opt){
  if(anal.type == "onedata"){
    SaveExpressClusterJSON(fileNm, clustOpt,opt);
  }else{
    initmetaloading <<- TRUE;
    SaveMetaClusterJSON(fileNm, clustOpt, opt);
  }
}

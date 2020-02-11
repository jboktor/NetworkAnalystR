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
#' @rdname GetExpressResultColNames
#' @export 
GetExpressResultColNames<-function(){
  resT <- readRDS("ExpressResT.rda");
  colnames(resT);
}

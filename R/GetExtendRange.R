#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vec PARAM_DESCRIPTION
#' @param unit PARAM_DESCRIPTION, Default: 10
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetExtendRange
#' @export 
GetExtendRange<-function(vec, unit=10){
  var.max <- max(vec);
  var.min <- min(vec);
  exts <- (var.max - var.min)/unit;
  c(var.min-exts, var.max+exts);
}

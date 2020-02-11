#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param boxplotName PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotDataBox
#' @export 
PlotDataBox <- function(boxplotName, dpi, format){
  qc.boxplot(dataSet$data.norm, boxplotName, dpi, format);
}

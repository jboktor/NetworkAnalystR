#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pcaName PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION
#' @param factor PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotDataPCA
#' @export 
PlotDataPCA <- function(pcaName, dpi, format,factor){
  qc.pcaplot(dataSet$data.norm, pcaName, dpi, format, factor);
}

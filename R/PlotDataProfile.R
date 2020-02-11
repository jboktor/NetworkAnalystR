#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param boxplotName PARAM_DESCRIPTION
#' @param pcaName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotDataProfile
#' @export 
PlotDataProfile<-function(dataName, boxplotName, pcaName){
  dataSet <- readRDS(dataName);
  qc.boxplot2(dataSet$data, boxplotName);
  qc.pcaplot2(dataSet$data, pcaName);
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dat PARAM_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname qc.meanstd
#' @export 
qc.meanstd <- function(dat, imgNm,dpi=72, format="png"){
  dpi = as.numeric(dpi)
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  plot = meanSdPlot(dat, ranks=FALSE) 
  dev.off();
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dat PARAM_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param format PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname qc.densitySample
#' @export 
qc.densitySample<- function(dat, imgNm, dpi=72, format){
  library("ggplot2")
  
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi = as.numeric(dpi)
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  df = data.frame(dataSet$data.norm, stringsAsFactors = FALSE)
  df = stack(df)
  g = ggplot(df, aes(x=values, color=ind)) +
    geom_density()
  
  print(g)
  dev.off();
}

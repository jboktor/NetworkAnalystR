#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
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
#' @rdname qc.metaDensity
#' @export 
qc.metaDensity<- function(imgNm, dpi=72, format, factor){
  library("ggplot2")
  inmex.meta <- readRDS("inmex_meta.rds");
  dat = inmex.meta$data
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi = as.numeric(dpi)
  
  df = data.frame(inmex.meta$data, stringsAsFactors = FALSE)
  df = stack(df)
  if(factor != "NA"){
    #Factor = as.vector(dataSet$meta.info[,factor])
  }else{
    #Factor = dataSet$meta.info[,1];
  }
  Factor=inmex.meta$data.lbl
  
  conv = data.frame(ind=colnames(inmex.meta$data), class=Factor)
  conv$ind=gsub("-", ".", conv$ind)
  df1 = merge(df, conv, by="ind")
  Cairo(file=imgNm, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  g =ggplot(df1, aes(x=values)) + geom_line(aes(color=class, group=ind), stat="density", alpha=0.3) + geom_line(aes(color=class), stat="density", alpha=0.6, size=1.5)
  print(g)
  dev.off();
}

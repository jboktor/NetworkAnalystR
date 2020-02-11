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
#' @rdname qc.density
#' @export 
qc.density<- function(imgNm, dpi=72, format, factor){
  library("ggplot2")
  dat = dataSet$data.norm
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi = as.numeric(dpi)
  
  df = data.frame(dataSet$data.norm, stringsAsFactors = FALSE)
  df = stack(df)
  sampleNms =gsub("-", ".", colnames(dataSet$data.norm))
  if(length(dataSet$meta.info) == 2){
    Factor1 = as.vector(dataSet$meta.info[,1])
    factorNm1 = colnames(dataSet$meta.info)[1]
    conv = data.frame(ind=sampleNms, Factor1=Factor1)
    colnames(conv) = c("ind", factorNm1);
    df1 = merge(df, conv, by="ind")
    Factor2 = as.vector(dataSet$meta.info[,2])
    factorNm2 = colnames(dataSet$meta.info)[2]
    conv = data.frame(ind=sampleNms, Factor2=Factor2)
    colnames(conv) = c("ind", factorNm2);
    df1 = merge(df1, conv, by="ind")
    df2 <- melt(df1, measure.vars=c(factorNm1,factorNm2))
    colnames(df2)[4] = "Conditions"
    g =ggplot(df2, aes(x=values)) + geom_line(aes(color=Conditions, group=ind), stat="density", alpha=0.6) + facet_grid(. ~ variable)
    width = 12
    height = 6
  }else{
    Conditions= as.character(dataSet$meta.info[,1]);
    conv = data.frame(ind=sampleNms, Conditions=Conditions)
    df1 = merge(df, conv, by="ind")
    g =ggplot(df1, aes(x=values)) + geom_line(aes(color=Conditions, group=ind), stat="density", alpha=0.6) 
    width = 8
    height = 6
  }
  Cairo(file=imgNm, width=width, height=height, type=format, bg="white", dpi=dpi, unit="in");
  print(g)
  dev.off();
}

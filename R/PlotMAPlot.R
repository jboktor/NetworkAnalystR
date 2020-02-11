#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param format PARAM_DESCRIPTION
#' @param pvalue PARAM_DESCRIPTION
#' @param fc PARAM_DESCRIPTION
#' @param inx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotMAPlot
#' @export 
PlotMAPlot<- function(imgNm, dpi=72, format,pvalue, fc, inx){
  library('ggplot2');
  dpi=as.numeric(dpi);
  inx = as.numeric(inx);
  pvalue=as.numeric(pvalue);
  fc=as.numeric(fc);
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=700, height=560,type=format, bg="white",dpi=100)
  pvalue = as.numeric(pvalue);
  fc = as.numeric(fc);
  res = dataSet$resTable
  if(dataSet$de.method =="deseq2"){
    #res = res[!apply(sapply(res, function(x) abs(scale(x)) >= 5), 1, any), ];
    res['log2(baseMean)'] = log2(res$baseMean);
  }
  
  # select based on p-value
  if(dataSet$type == "array"){
    res$significant = ifelse(res[,inx] > fc & res$P.Value < pvalue , 1, ifelse(res[,inx] < -fc & res$P.Value < pvalue, -1, 0))
  } else {
    res$significant = ifelse(res[,inx] > fc & res$adj.P.Val < pvalue , 1, ifelse(res[,inx] < -fc & res$adj.P.Val < pvalue, -1, 0))
  }
  res$significant = as.factor(res$significant)
  yCol= colnames(res)[inx]
  
  if (dataSet$de.method=="limma"){
    maplot = ggplot(res, aes_string(x="AveExpr", y=yCol, color="significant"))
  } else if (dataSet$de.method=="deseq2"){
    maplot = ggplot(res, aes_string(x="log2(baseMean)", y=yCol, color="significant"))
  } else {
    maplot = ggplot(res, aes_string(x="logCPM", y=yCol, color="significant"))
  }
  
  maplot = maplot +
    geom_point(size=0.5, alpha=0.5) +
    geom_hline(color = "blue3", yintercept = 0) +
    stat_smooth(se = FALSE, method = "loess", color = "red3") +
    scale_color_manual(values=c("-1" = "green", "0" = "black", "1" = "red"))+ theme(legend.position="none")
  print(maplot)
  dev.off()
}

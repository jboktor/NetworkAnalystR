#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
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
#' @rdname PlotMetaPCA
#' @export 
PlotMetaPCA <- function(imgNm, dpi, format,factor){
  inmex.meta <- readRDS("inmex_meta.rds");
  x=inmex.meta[["data"]]
  dpi = as.numeric(dpi);
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  library('lattice');
  library('ggplot2');
  pca <- prcomp(t(na.omit(x)));
  imp.pca<-summary(pca)$importance;
  xlabel = paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
  ylabel = paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  if(factor != "NA"){
    #Factor = as.vector(dataSet$meta.info[,factor])
  }else{
    #Factor = dataSet$meta.info[,1];
  }
  Conditions = factor(inmex.meta$cls.lbl)
  Datasets = factor(inmex.meta$data.lbl)
  pcafig = ggplot(pca.res, aes(x=PC1, y=PC2,  color=Conditions ,shape=Datasets)) +
    geom_point(size=4, alpha=0.5) + xlim(xlim)+ ylim(ylim) + xlab(xlabel) + ylab(ylabel);
  
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", unit="in", dpi=dpi);
  print(pcafig);
  dev.off();
  
}

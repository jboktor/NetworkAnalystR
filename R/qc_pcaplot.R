#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @param factor PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname qc.pcaplot
#' @export 
qc.pcaplot <- function(x, imgNm, dpi=72, format="png", factor){
  dpi = as.numeric(dpi);
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  require('lattice');
  require('ggplot2');
  require('reshape');
  pca <- prcomp(t(na.omit(x)));
  imp.pca<-summary(pca)$importance;
  xlabel = paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
  ylabel = paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  pca.res <- pca.res[,c(1,2)]
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  if(length(dataSet$meta.info) == 2){
    Factor1 = as.vector(dataSet$meta.info[,1])
    factorNm1 = colnames(dataSet$meta.info)[1]
    pca.res[,factorNm1] = Factor1
    Factor2 = as.vector(dataSet$meta.info[,2])
    factorNm2 = colnames(dataSet$meta.info)[2]
    pca.res[,factorNm2] = Factor2
    pca.rest <- melt(pca.res, measure.vars=c(factorNm1,factorNm2))
    colnames(pca.rest)[4] = "Conditions"
    pca.rest$names = c(rownames(pca.res), rownames(pca.res))
    if(length(pca.rest$names)>20){
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions, label=pca.rest$names)) +
        geom_point(size=3, alpha=0.5) + xlim(xlim)+ ylim(ylim) + xlab(xlabel) + ylab(ylabel) + facet_grid(. ~ variable)
    }else{
      require('ggrepel');
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions, label=pca.rest$names)) +
        geom_point(size=4) + xlim(xlim)+ ylim(ylim) + xlab(xlabel) + ylab(ylabel) + geom_text_repel(force=1.5) + facet_grid(. ~ variable)
    }
    width = 12
    height = 6
  }else{
    Factor = dataSet$meta.info[,1];
    pca.rest = pca.res
    pca.rest$Conditions = Factor
    pca.rest$names = rownames(pca.res)
    if(length(rownames(pca.res))>20){
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
        geom_point(size=3, alpha=0.5) + xlim(xlim)+ ylim(ylim)+ xlab(xlabel) + ylab(ylabel) 
    }else{
      require('ggrepel');
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions, label=rownames(pca.res))) +
        geom_point(size=4) + xlim(xlim)+ ylim(ylim)+ xlab(xlabel) + ylab(ylabel) +geom_text_repel(force=1.5)+scale_color_manual(breaks=unique(pca.rest$Conditions), values=c("#00BFC4" ,"#F8766D"))
    }
    width = 8
    height = 6
  }
  Cairo(file=imgNm, width=width, height=height, type=format, bg="white", unit="in", dpi=dpi);
  print(pcafig);
  dev.off();
}

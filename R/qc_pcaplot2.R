#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname qc.pcaplot2
#' @export 
qc.pcaplot2 <- function(x, imgNm){
  imgNm = paste(imgNm, "dpi", "72", ".png", sep="");
  require('lattice');
  pca <- prcomp(t(na.omit(x)));
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  pcafig = xyplot(PC2~PC1, data=pca.res, pch=19, cex=1,aspect = "iso", xlim = xlim, ylim=ylim,
                  panel=function(x, y, ...) {
                    panel.xyplot(x, y, ...);
                    ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.8, col="magenta")
                  })
  
  Cairo(file=imgNm, width=480, height=480, type="png", bg="white");
  print(pcafig);
  dev.off();
}

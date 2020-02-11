#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imgName PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotMDS
#' @export 
PlotMDS <- function(imgName, format){
  require(edgeR)
  require(RColorBrewer)
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgName, width=580, type=format, bg="white",dpi=72)
  levels(dataSet$cls) <- brewer.pal(nlevels(dataSet$cls), "Set1")
  col.group <- dataSet$cls
  col.group <- as.character(col.group)
  plotMDS(dataSet$data.norm, col=col.group, xlab = "Dimension 2", ylab = "Dimension 1")
  title(main="MDS")
  dev.off();
} 

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cmpdNm PARAM_DESCRIPTION
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param width PARAM_DESCRIPTION, Default: NA
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotGSView
#' @export 
PlotGSView <-function(cmpdNm,  format="png", dpi=72, width=NA){
  library("ggplot2");
  current.geneset <- readRDS("current_geneset.rds")
  nm <<- cmpdNm
  imgName <- gsub("\\/", "_",  cmpdNm);
  imgName <- gsub(" ", "_",  imgName);
  imgName <- paste(imgName, "_dpi", dpi, ".", format, sep="");
  #indx<-which(rownames(boxplot_id)==cmpdNm);
  #gene.id <- boxplot_id[indx,1];
  Cairo(file = imgName, dpi=72, width=340, height=300, type="png", bg="transparent");
  g <- plotEnrichment(current.geneset[[cmpdNm]], rankedVec)
  print(g)
  dev.off();
  return(imgName);
}

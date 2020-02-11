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
#' @rdname PlotCmpdView
#' @export 
PlotCmpdView <-function(cmpdNm, format="png", dpi=72, width=NA){
  if(anal.type == "onedata"){
    datanorm = dataSet$data.norm
  }else{
    datanorm = dataSet$data
  }
  clslbl = dataSet$meta.info[,1];
  imgName <- gsub("\\/", "_",  cmpdNm);
  imgName <- paste(imgName, "_dpi", dpi, ".", format, sep="");
  #indx<-which(rownames(boxplot_id)==cmpdNm);
  #gene.id <- boxplot_id[indx,1];
  gene.symb <<- doEntrez2SymbolMapping(cmpdNm);
  Cairo(file = imgName, dpi=dpi, width=230, height=230, type=format, bg="transparent");
  par(mar=c(4,3,1,2), oma=c(0,0,1,0));
  boxplot(datanorm[which(rownames(datanorm)==as.character(cmpdNm)),]~clslbl,las=2,col= unique(GetColorSchema(clslbl)));
  title(main=gene.symb, out=T);
  dev.off();
  return(imgName);
}

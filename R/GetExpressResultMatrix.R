#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param inxt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetExpressResultMatrix
#' @export 
GetExpressResultMatrix <-function(inxt){
  
  inxt = as.numeric(inxt)
  if (dataSet$de.method=="limma"){
    inx =match("AveExpr", colnames(dataSet$resTable))
  } else if (dataSet$de.method=="deseq2"){
    inx =match("baseMean", colnames(dataSet$resTable))
  } else {
    inx =match("logCPM", colnames(dataSet$resTable))
  }
  res = dataSet$resTable;
  res = res[,-(1:inx-1)]
  res <- cbind(dataSet$resTable[,inxt], res);
  colnames(res)[1] <- colnames(dataSet$resTable)[inxt];
  saveRDS(res, "ExpressResT.rda");
  return(signif(as.matrix(res), 5));
}

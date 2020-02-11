#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMetaCol
#' @export 
GetMetaCol<- function(){
  colNms <- colnames(dataSet$resTable);
  if (dataSet$de.method=="limma"){
    inx =match("AveExpr", colNms)
  } else if (dataSet$de.method=="deseq2"){
    inx =match("baseMean", colNms)
  } else {
    inx =match("logCPM", colNms)
  }
  resT <- dataSet$resTable;
  resT <- resT[,1:inx-1]
  return(colnames(resT));
}

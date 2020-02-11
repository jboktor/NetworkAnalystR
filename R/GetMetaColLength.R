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
#' @rdname GetMetaColLength
#' @export 
GetMetaColLength<- function(){
  if (dataSet$de.method=="limma"){
    inx =match("AveExpr", colnames(dataSet$resTable))
  } else if (dataSet$de.method=="deseq2"){
    inx =match("baseMean", colnames(dataSet$resTable))
  } else {
    inx =match("logCPM", colnames(dataSet$resTable))
  }
  resT = dataSet$resTable;
  resT = resT[,1:inx-1]
  return(length(colnames(resT)));
}

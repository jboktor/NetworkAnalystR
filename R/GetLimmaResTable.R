# get result table from eBayes fit object
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fit.obj PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetLimmaResTable
#' @export 
GetLimmaResTable<-function(fit.obj){
  resTable <- topTable(fit.obj, number=Inf, adjust.method="BH");
  if(!is.null(resTable$ID)){ # for older version
    rownames(resTable) <- resTable$ID;
    resTable$ID <- NULL;
  }
  return (resTable);
}

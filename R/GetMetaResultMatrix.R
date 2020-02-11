# single.type return logFC or p value for individual data analysis
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param single.type PARAM_DESCRIPTION, Default: 'fc'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetMetaResultMatrix
#' @export 
GetMetaResultMatrix<-function(single.type="fc"){
  if(single.type == "fc"){
    meta.mat <- cbind(fc.mat, meta.mat);
  }else{
    meta.mat <- cbind(pval.mat, meta.mat);
  }
  # display at most 500 genes
  if(nrow(meta.mat) > 500){
    meta.mat <- meta.mat[1:500,]; # already sorted based on meta-p values
  }
  meta.mat <-signif(as.matrix(meta.mat), 5);
  meta.mat;
}

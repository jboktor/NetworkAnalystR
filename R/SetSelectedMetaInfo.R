#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param meta0 PARAM_DESCRIPTION
#' @param meta1 PARAM_DESCRIPTION
#' @param block1 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetSelectedMetaInfo
#' @export 
SetSelectedMetaInfo <- function(meta0, meta1, block1){
  if(meta0 == "NA"){
    return(0);
  }else{
    cls <- dataSet$meta.info[, meta0];
    dataSet$fst.cls <- cls; # for PCA plotting
    block <- NULL;
    dataSet$sec.cls <- "NA"
    if(meta1 != "NA"){
      if(block1){
        block <- dataSet$meta.info[, meta1]
      }else{ # two factor
        cls <- interaction(dataSet$meta.info[, c(meta0, meta1)], sep = ".", lex.order = TRUE);
      }
      dataSet$sec.cls <- dataSet$meta.info[, meta1]; # for pca coloring
    }
    dataSet$cls <- cls; # record main cls;
    dataSet$block <- block;
    dataSet <<- dataSet;
    return(levels(cls));
  }
}

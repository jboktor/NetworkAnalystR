# obtain sample names and their class labels
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param clsLbl PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetSampleInfo
#' @export 
GetSampleInfo <- function(dataName, clsLbl){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  grpInfo <- dataSet$meta.info[[clsLbl]];
  grpLbls <- paste(levels(grpInfo), collapse="\n");
  smplInfo <- paste(Sample = colnames(dataSet$data.orig), "\t", Class=grpInfo, collapse="\n");
  return(c(grpLbls, smplInfo));
}

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
#' @rdname GetMetaResultGeneIDs
#' @export 
GetMetaResultGeneIDs<-function(){
  rnms <- rownames(as.matrix(meta.mat));# already sorted based on meta-p values
  if(length(rnms) > 500){
    rnms <- rnms[1:500];
  }
  return(rnms);
}

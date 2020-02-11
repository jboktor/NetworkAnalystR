#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ids PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetVolcanoHigh
#' @export 
SetVolcanoHigh<- function(ids){
  idsu <<- ids
  gene.vec <- unlist(strsplit(ids, "; "));
  volcanoHlVec <<- gene.vec;
  if(length(volcanoHlVec)>0){
    return(1);
  }else{
    return(0);
  }
}

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
#' @rdname GetQEA.pathNames
#' @export 
GetQEA.pathNames<-function(){
  current.geneset <- readRDS("current_geneset.rds")
  hit.inx <- match(rownames(analSet$qea.mat),names(current.geneset));
  return(names(current.geneset)[hit.inx]);
}

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
#' @rdname GetMetaGeneIDType
#' @export 
GetMetaGeneIDType<-function(){
  inmex.meta <- readRDS("inmex_meta.rds");
  return(inmex.meta$id.type);
}

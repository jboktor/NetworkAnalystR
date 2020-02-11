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
#' @rdname GetSelectedDataGeneNumberUpdated
#' @export 
GetSelectedDataGeneNumberUpdated<- function(){
  return(paste(venn.genenb.up, collapse=";"));
}

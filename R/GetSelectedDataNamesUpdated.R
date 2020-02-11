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
#' @rdname GetSelectedDataNamesUpdated
#' @export 
GetSelectedDataNamesUpdated <- function(){
  return(paste(names(venn.list.up), collapse=";"));
}

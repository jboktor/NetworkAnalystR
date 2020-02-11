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
#' @rdname GetSelectedDataGeneNumber
#' @export 
GetSelectedDataGeneNumber<- function(){
  return(paste(venn.genenb, collapse=";"));
}

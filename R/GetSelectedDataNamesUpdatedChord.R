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
#' @rdname GetSelectedDataNamesUpdatedChord
#' @export 
GetSelectedDataNamesUpdatedChord <- function(){
  return(paste(names(chord.list.up), collapse=";"));
}

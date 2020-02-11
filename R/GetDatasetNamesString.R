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
#' @rdname GetDatasetNamesString
#' @export 
GetDatasetNamesString <- function(){
  inmex.meta <- readRDS("inmex_meta.rds");
  paste(unique(inmex.meta$data.lbl), collapse="||");
}
